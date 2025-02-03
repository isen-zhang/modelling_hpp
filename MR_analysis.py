import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from typing import Tuple, List, Dict
import matplotlib.pyplot as plt
import seaborn as sns
import os

def load_gwas_data(phenotype: str, analysis_type: str) -> pd.DataFrame:
    """
    Load GWAS results for a specific phenotype
    
    Parameters:
    -----------
    phenotype: str
        Name of the phenotype
    analysis_type: str
        Type of analysis ('linear' or 'logistic')
    
    Returns:
    --------
    pd.DataFrame
        Processed GWAS results
    """
    # Construct file path
    filepath = f"gwas_results/gwas_{phenotype}.assoc.{analysis_type}"
    
    # Read the GWAS results
    gwas_data = pd.read_csv(filepath, delim_whitespace=True)
    
    # Print debug info
    print(f"Loading {analysis_type} GWAS data for {phenotype}")
    print("Original columns:", gwas_data.columns.tolist())
    
    # Convert OR to BETA if it's logistic regression
    if analysis_type == 'logistic':
        if 'OR' in gwas_data.columns:
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
    
    # Rename only the columns that exist
    existing_columns = {k: v for k, v in column_mapping.items() if k in gwas_data.columns}
    gwas_data = gwas_data.rename(columns=existing_columns)
    
    # Calculate SE based on the analysis type
    if analysis_type == 'logistic':
        # For logistic regression, STAT is typically Z-score
        gwas_data['se'] = abs(gwas_data['beta'] / gwas_data['stat'])
    else:
        # For linear regression
        gwas_data['se'] = abs(gwas_data['beta'] / gwas_data['stat'])
    
    return gwas_data

class GWASMendelianRandomization:
    def __init__(self, 
                 exposure_phenotype: str,
                 outcome_phenotype: str,
                 exposure_type: str = 'linear',
                 outcome_type: str = 'linear',
                 pval_threshold: float = 5e-8):
        """
        Initialize MR analysis with GWAS results
        
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
        print(f"\nLoading exposure data: {exposure_phenotype} ({exposure_type})")
        self.exposure_data = load_gwas_data(exposure_phenotype, exposure_type)
        
        print(f"\nLoading outcome data: {outcome_phenotype} ({outcome_type})")
        self.outcome_data = load_gwas_data(outcome_phenotype, outcome_type)
        
        # Select instrumental variables
        self.exposure_data = self.exposure_data[
            self.exposure_data['pval'] < pval_threshold
        ].copy()
        
        print(f"\nSelected {len(self.exposure_data)} SNPs as instruments (p < {pval_threshold})")
        
        # Store analysis types
        self.exposure_type = exposure_type
        self.outcome_type = outcome_type
        
        # Initialize harmonized data
        self.harmonized_data = None
    
    def harmonize_data(self) -> pd.DataFrame:
        """Harmonize the exposure and outcome data"""
        # Merge datasets
        merged_data = pd.merge(
            self.exposure_data,
            self.outcome_data,
            on='SNP',
            suffixes=('_exposure', '_outcome')
        )
        
        print(f"\nHarmonized {len(merged_data)} SNPs between exposure and outcome")
        
        # Ensure effect alleles are aligned
        allele_flip = merged_data['effect_allele_exposure'] != merged_data['effect_allele_outcome']
        merged_data.loc[allele_flip, 'beta_outcome'] *= -1
        
        print(f"Flipped effect alleles for {sum(allele_flip)} SNPs")
        
        # Store harmonized data
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
    
    def create_diagnostic_plots(self, output_dir: str = 'mr_plots'):
        """Create diagnostic plots for MR analysis"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Scatter plot
        plt.figure(figsize=(10, 8))
        plt.scatter(self.harmonized_data['beta_exposure'],
                   self.harmonized_data['beta_outcome'],
                   alpha=0.6)
        plt.xlabel(f'SNP effect on exposure ({self.exposure_type})')
        plt.ylabel(f'SNP effect on outcome ({self.outcome_type})')
        plt.title('MR Scatter Plot')
        
        # Add IVW slope
        beta_ivw, _, _ = self.ivw_method()
        x_range = np.array([
            self.harmonized_data['beta_exposure'].min(),
            self.harmonized_data['beta_exposure'].max()
        ])
        plt.plot(x_range, beta_ivw * x_range, 'r--', label=f'IVW (β={beta_ivw:.3f})')
        plt.legend()
        plt.savefig(f'{output_dir}/scatter_plot.png')
        plt.close()
        
        # Funnel plot
        plt.figure(figsize=(10, 8))
        ratios = self.harmonized_data['beta_outcome'] / self.harmonized_data['beta_exposure']
        precision = 1 / self.harmonized_data['se_outcome']
        plt.scatter(ratios, precision, alpha=0.6)
        plt.xlabel('SNP effect estimate')
        plt.ylabel('Precision')
        plt.title('Funnel Plot')
        plt.savefig(f'{output_dir}/funnel_plot.png')
        plt.close()
    
    def run_analysis(self) -> Dict:
        """Run complete MR analysis"""
        # Harmonize data
        self.harmonize_data()
        
        # Run IVW method
        ivw_beta, ivw_se, ivw_pval = self.ivw_method()
        
        # Create diagnostic plots
        self.create_diagnostic_plots()
        
        # Compile results
        results = {
            'IVW': {
                'beta': ivw_beta,
                'se': ivw_se,
                'pval': ivw_pval,
                'n_snps': len(self.harmonized_data)
            },
            'analysis_info': {
                'exposure_type': self.exposure_type,
                'outcome_type': self.outcome_type,
                'n_exposure_snps': len(self.exposure_data),
                'n_harmonized_snps': len(self.harmonized_data)
            }
        }
        
        return results

# Example usage:
"""
# Initialize MR analysis
mr_analysis = GWASMendelianRandomization(
    exposure_phenotype='height',
    outcome_phenotype='bmi',
    exposure_type='logistic',  # or 'linear'
    outcome_type='linear',     # or 'logistic'
    pval_threshold=5e-8
)

# Run analysis
results = mr_analysis.run_analysis()

# Print results
print("\nMR Analysis Results:")
print(f"Number of exposure SNPs: {results['analysis_info']['n_exposure_snps']}")
print(f"Number of harmonized SNPs: {results['analysis_info']['n_harmonized_snps']}")
print("\nIVW Results:")
print(f"Beta: {results['IVW']['beta']:.3f}")
print(f"SE: {results['IVW']['se']:.3f}")
print(f"P-value: {results['IVW']['pval']:.3e}")
print(f"Number of SNPs used: {results['IVW']['n_snps']}")
"""
