import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from typing import Tuple, Dict
import matplotlib.pyplot as plt
import seaborn as sns
import os

class ComprehensiveMR:
    def __init__(self, 
                 exposure_phenotype: str,
                 outcome_phenotype: str,
                 exposure_type: str = 'linear',
                 outcome_type: str = 'linear',
                 pval_threshold: float = 5e-8):
        """
        Initialize comprehensive MR analysis
        
        Parameters:
        -----------
        exposure_phenotype: str
            Name of the exposure phenotype
        outcome_phenotype: str
            Name of the outcome phenotype
        exposure_type: str
            Type of GWAS analysis for exposure ('linear' or 'logistic')
        outcome_type: str
            Type of GWAS analysis for outcome ('linear' or 'logistic')
        pval_threshold: float
            P-value threshold for selecting instrumental variables
        """
        self.exposure_phenotype = exposure_phenotype
        self.outcome_phenotype = outcome_phenotype
        self.exposure_type = exposure_type
        self.outcome_type = outcome_type
        
        print(f"\nLoading exposure data: {exposure_phenotype} ({exposure_type})")
        self.exposure_data = self.load_gwas_data(exposure_phenotype, exposure_type)
        
        print(f"\nLoading outcome data: {outcome_phenotype} ({outcome_type})")
        self.outcome_data = self.load_gwas_data(outcome_phenotype, outcome_type)
        
        # Select instrumental variables
        self.exposure_data = self.exposure_data[
            self.exposure_data['pval'] < pval_threshold
        ].copy()
        
        print(f"\nSelected {len(self.exposure_data)} SNPs as instruments (p < {pval_threshold})")
        
        # Initialize harmonized data
        self.harmonized_data = None
    
    @staticmethod
    def load_gwas_data(phenotype: str, analysis_type: str) -> pd.DataFrame:
        """Load and process GWAS results"""
        filepath = f"gwas_results/gwas_{phenotype}.assoc.{analysis_type}"
        gwas_data = pd.read_csv(filepath, delim_whitespace=True)
        
        print("Original columns:", gwas_data.columns.tolist())
        
        # Convert OR to BETA for logistic regression
        if analysis_type == 'logistic' and 'OR' in gwas_data.columns:
            gwas_data['BETA'] = np.log(gwas_data['OR'])
            print("Converted OR to BETA using natural logarithm")
        
        # Standardize column names
        column_mapping = {
            'SNP': 'SNP',
            'BP': 'BP',
            'A1': 'effect_allele',
            'TEST': 'TEST',
            'NMISS': 'N',
            'BETA': 'beta',
            'STAT': 'stat',
            'P': 'pval'
        }
        
        # Rename only existing columns
        existing_columns = {k: v for k, v in column_mapping.items() if k in gwas_data.columns}
        gwas_data = gwas_data.rename(columns=existing_columns)
        
        # Calculate SE
        gwas_data['se'] = abs(gwas_data['beta'] / gwas_data['stat'])
        
        return gwas_data
    
    def harmonize_data(self) -> pd.DataFrame:
        """Harmonize the exposure and outcome data"""
        merged_data = pd.merge(
            self.exposure_data,
            self.outcome_data,
            on='SNP',
            suffixes=('_exposure', '_outcome')
        )
        
        print(f"\nHarmonized {len(merged_data)} SNPs between exposure and outcome")
        
        allele_flip = merged_data['effect_allele_exposure'] != merged_data['effect_allele_outcome']
        merged_data.loc[allele_flip, 'beta_outcome'] *= -1
        
        print(f"Flipped effect alleles for {sum(allele_flip)} SNPs")
        
        self.harmonized_data = merged_data
        return self.harmonized_data
    
    def ivw_method(self) -> Tuple[float, float, float]:
        """Perform inverse variance weighted method"""
        weights = 1 / (self.harmonized_data['se_outcome'] ** 2)
        
        beta = np.sum(weights * self.harmonized_data['beta_outcome'] * 
                     self.harmonized_data['beta_exposure']) / \
               np.sum(weights * self.harmonized_data['beta_exposure'] ** 2)
        
        se = np.sqrt(1 / np.sum(weights * self.harmonized_data['beta_exposure'] ** 2))
        
        z_score = beta / se
        pval = 2 * (1 - stats.norm.cdf(abs(z_score)))
        
        return beta, se, pval
    
    def mr_egger(self) -> Dict:
        """Perform MR-Egger regression"""
        X = sm.add_constant(self.harmonized_data['beta_exposure'])
        y = self.harmonized_data['beta_outcome']
        weights = 1 / (self.harmonized_data['se_outcome'] ** 2)
        
        model = sm.WLS(y, X, weights=weights)
        results = model.fit()
        
        return {
            'slope': results.params[1],
            'slope_se': results.bse[1],
            'slope_pval': results.pvalues[1],
            'intercept': results.params[0],
            'intercept_se': results.bse[0],
            'intercept_pval': results.pvalues[0]
        }
    
    def weighted_median(self) -> Dict:
        """Perform weighted median estimation"""
        weights = 1 / (self.harmonized_data['se_outcome'] ** 2)
        beta_iv = self.harmonized_data['beta_outcome'] / self.harmonized_data['beta_exposure']
        
        weights = weights / np.sum(weights)
        
        sorted_indices = np.argsort(beta_iv)
        sorted_weights = weights[sorted_indices]
        sorted_betas = beta_iv[sorted_indices]
        
        cumulative_weights = np.cumsum(sorted_weights)
        median_index = np.searchsorted(cumulative_weights, 0.5)
        estimate = sorted_betas[median_index]
        
        # Bootstrap for SE
        n_boot = 1000
        boot_estimates = []
        for _ in range(n_boot):
            boot_indices = np.random.choice(
                len(beta_iv),
                size=len(beta_iv),
                p=weights
            )
            boot_estimates.append(np.median(beta_iv[boot_indices]))
        
        se = np.std(boot_estimates)
        pval = 2 * (1 - stats.norm.cdf(abs(estimate / se)))
        
        return {
            'beta': estimate,
            'se': se,
            'pval': pval
        }
    
    def leave_one_out(self) -> pd.DataFrame:
        """Perform leave-one-out analysis"""
        results = []
        for i in range(len(self.harmonized_data)):
            subset = self.harmonized_data.drop(self.harmonized_data.index[i])
            weights = 1 / (subset['se_outcome'] ** 2)
            beta = (np.sum(weights * subset['beta_outcome'] * subset['beta_exposure']) /
                   np.sum(weights * subset['beta_exposure'] ** 2))
            se = np.sqrt(1 / np.sum(weights * subset['beta_exposure'] ** 2))
            
            results.append({
                'excluded_snp': self.harmonized_data.iloc[i]['SNP'],
                'beta': beta,
                'se': se,
                'pval': 2 * (1 - stats.norm.cdf(abs(beta / se)))
            })
        
        return pd.DataFrame(results)
    
    def create_plots(self, output_dir: str = 'mr_diagnostic_plots'):
        """Create diagnostic plots"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Format phenotype names for titles
        exposure_title = self.exposure_phenotype.replace('_', ' ').title()
        outcome_title = self.outcome_phenotype.replace('_', ' ').title()
        base_filename = f"{self.exposure_phenotype}_to_{self.outcome_phenotype}"
        
        # Scatter plot
        plt.figure(figsize=(10, 8))
        plt.scatter(self.harmonized_data['beta_exposure'], 
                   self.harmonized_data['beta_outcome'],
                   alpha=0.6, label='SNPs')
        
        # Add regression lines
        x_range = np.array([
            self.harmonized_data['beta_exposure'].min(),
            self.harmonized_data['beta_exposure'].max()
        ])
        
        # IVW line
        ivw_beta, _, _ = self.ivw_method()
        plt.plot(x_range, ivw_beta * x_range, 'r--',
                label=f'IVW (β={ivw_beta:.3f})')
        
        # MR-Egger line
        egger_results = self.mr_egger()
        plt.plot(x_range,
                egger_results['intercept'] + egger_results['slope'] * x_range,
                'g--', label=f'MR-Egger (β={egger_results["slope"]:.3f})')
        
        plt.xlabel(f'SNP effect on {exposure_title}')
        plt.ylabel(f'SNP effect on {outcome_title}')
        plt.title(f'MR Scatter Plot: {exposure_title} → {outcome_title}')
        plt.legend()
        plt.savefig(f'{output_dir}/{base_filename}_scatter_plot.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Funnel plot
        plt.figure(figsize=(10, 8))
        precision = 1 / self.harmonized_data['se_outcome']
        beta_iv = self.harmonized_data['beta_outcome'] / self.harmonized_data['beta_exposure']
        plt.scatter(beta_iv, precision, alpha=0.6)
        plt.axvline(x=ivw_beta, color='r', linestyle='--',
                   label=f'IVW estimate: {ivw_beta:.3f}')
        plt.xlabel('SNP-specific causal estimate')
        plt.ylabel('Precision (1/SE)')
        plt.title(f'Funnel Plot: {exposure_title} → {outcome_title}')
        plt.legend()
        plt.savefig(f'{output_dir}/{base_filename}_funnel_plot.png',
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Leave-one-out plot
        loo_results = self.leave_one_out()
        
        # Calculate figure height based on number of SNPs
        height_per_snp = 0.25  # inches per SNP
        fig_height = max(6, len(loo_results) * height_per_snp)
        
        plt.figure(figsize=(12, fig_height))
        
        # Create the error bar plot
        plt.errorbar(loo_results['beta'], range(len(loo_results)),
                    xerr=1.96 * loo_results['se'], fmt='o', alpha=0.5)
        
        # Add SNP names on y-axis
        plt.yticks(range(len(loo_results)), loo_results['excluded_snp'], fontsize=8)
        
        # Add IVW estimate line
        plt.axvline(x=ivw_beta, color='r', linestyle='--',
                   label=f'IVW estimate (β={ivw_beta:.3f})')
        
        plt.xlabel('Causal estimate')
        plt.title(f'Leave-one-out Analysis: {exposure_title} → {outcome_title}')
        plt.legend()
        
        # Adjust layout to prevent SNP names from being cut off
        plt.tight_layout()
        
        # Save with larger height if needed
        plt.savefig(f'{output_dir}/{base_filename}_leave_one_out_plot.png',
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    def run_analysis(self) -> Dict:
        """Run complete MR analysis"""
        # Harmonize data if not already done
        if self.harmonized_data is None:
            self.harmonize_data()
        
        # Run all methods
        ivw_beta, ivw_se, ivw_pval = self.ivw_method()
        egger_results = self.mr_egger()
        weighted_median_results = self.weighted_median()
        
        # Create diagnostic plots
        self.create_plots()
        
        # Compile results
        results = {
            'IVW': {
                'beta': ivw_beta,
                'se': ivw_se,
                'pval': ivw_pval,
                'n_snps': len(self.harmonized_data)
            },
            'MR_Egger': egger_results,
            'Weighted_Median': weighted_median_results,
            'analysis_info': {
                'exposure_phenotype': self.exposure_phenotype,
                'outcome_phenotype': self.outcome_phenotype,
                'exposure_type': self.exposure_type,
                'outcome_type': self.outcome_type,
                'n_exposure_snps': len(self.exposure_data),
                'n_harmonized_snps': len(self.harmonized_data)
            }
        }
        
        return results

# Example usage:
"""
# Initialize and run analysis
mr = ComprehensiveMR(
    exposure_phenotype='height',
    outcome_phenotype='back_pain',
    exposure_type='linear',
    outcome_type='logistic',
    pval_threshold=5e-8
)

# Run analysis
results = mr.run_analysis()

# Print results
print("\nMR Analysis Results:")
print(f"Exposure: {results['analysis_info']['exposure_phenotype']}")
print(f"Outcome: {results['analysis_info']['outcome_phenotype']}")
print(f"\nNumber of exposure SNPs: {results['analysis_info']['n_exposure_snps']}")
print(f"Number of harmonized SNPs: {results['analysis_info']['n_harmonized_snps']}")

print("\nIVW Results:")
print(f"Beta: {results['IVW']['beta']:.3f}")
print(f"SE: {results['IVW']['se']:.3f}")
print(f"P-value: {results['IVW']['pval']:.3e}")

print("\nMR-Egger Results:")
print(f"Slope: {results['MR_Egger']['slope']:.3f}")
print(f"Slope P-value: {results['MR_Egger']['slope_pval']:.3e}")
print(f"Intercept: {results['MR_Egger']['intercept']:.3f}")
print(f"Intercept P-value: {results['MR_Egger']['intercept_pval']:.3e}")

print("\nWeighted Median Results:")
print(f"Beta: {results['Weighted_Median']['beta']:.3f}")
print(f"P-value: {results['Weighted_Median']['pval']:.3e}")
"""
