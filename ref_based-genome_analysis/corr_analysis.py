import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import logging
from pathlib import Path
import numpy as np 

input_files = {
    'mincov10_sig': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score1_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score1_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score1_sum_v2+.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score2_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score2_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score2_sum_v2+.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score3_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score3_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score3_sum_v2+.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score4_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score4_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_0.05/score4_sum_v2+.csv'
        }
    },
    'mincov10_minedit0.05': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score_count_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score_sum_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score_mean_v1.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score2_count_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score2_sum_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score2_mean_v1.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score3_count_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score3_sum_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score3_mean_v1.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score4_count_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score4_sum_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/19apr_mincov10_alledits/sample_metrics_mincov10_minedit.05_dbSNP/v3_unfiltered/score4_mean_v1.csv'
        }
    },
    'mincov10_minedit0.05_sig': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score1_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score1_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score1_sum_v2+.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score2_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score2_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score2_sum_v2+.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score3_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score3_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score3_sum_v2+.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score4_count_v2+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score4_mean_v2+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_mincov10_minedit0.05_0.05/score4_sum_v2+.csv'
        }
    },
    'mincov1': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 1/rep_count_score1_0.01+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 1/rep_mean_score1_0.01+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 1/rep_sum_score1_0.01+.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 2/rep_count_score2_log+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 2/rep_mean_score2_log+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 2/rep_sum_score2_log+.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 3/rep_count_score3_zscore.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 3/rep_mean_score3_zscore.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 3/rep_sum_score3_zscore.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 4/rep_count_score4_zscore1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 4/rep_mean_score4_zscore1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/MINCOV1/04.sample_metrics_mincov1/score 4/rep_sum_score4_zscore1.csv'
            }
    },
    'mincov10': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score1/rep_count_score1_0.01+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score1/rep_mean_score1_0.01+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score1/rep_sum_score1_0.01+.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score2/rep_count_score2_log+.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score2/rep_mean_score2_log+.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score2/rep_sum_score2_log+.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score3/rep_count_score3_zscore.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score3/rep_mean_score3_zscore.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score3/rep_sum_score3_zscore.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score4/rep_count_score4_zscore1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score4/rep_mean_score4_zscore1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/bin/04.sample_metrics_final_v3/score4/rep_sum_score4_zscore1.csv'
        }
    },
    'nomincov_sig': {
        'score1': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score1_count_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score1_mean_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score1_sum_v1.csv'
        },
        'score2': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score2_count_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score2_mean_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score2_sum_v1.csv'
        },
        'score3': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score3_count_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score3_mean_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score3_sum_v1.csv'
        },
        'score4': {
            'count': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score4_count_v1.csv',
            'mean': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score4_mean_v1.csv',
            'sum': '/Volumes/Lucky me/8_STAMP/25Mar_fisher_0.05/score4_sum_v1.csv'
        }
    }
}

# Score column mappings
score_column_mappings = {
    'score1': {
        'primary': {
            'SKOV3_Apo_RIPK1_BT1_UT_ratio': 'ratio_RIPK1',
            'SKOV3_Apo_FTO_BT1_UT_ratio': 'ratio_ApoFTO', 
            'SKOV3_Tad_FTO_BT1_UT_ratio': 'ratio_TadFTO'
        },
        'alternative': {
            'ratio_RIPK1': 'ratio_RIPK1',
            'ratio_ApoFTO': 'ratio_ApoFTO',
            'ratio_TadFTO': 'ratio_TadFTO'
        }
    },
    'score2': {
        'primary': {
            'SKOV3_Apo_RIPK1_BT1_UT_ratio_log': 'ratio_log_RIPK1',
            'SKOV3_Apo_FTO_BT1_UT_ratio_log': 'ratio_log_ApoFTO', 
            'SKOV3_Tad_FTO_BT1_UT_ratio_log': 'ratio_log_TadFTO'
        },
        'alternative': {
            'ratio_RIPK1': 'ratio_log_RIPK1',
            'ratio_ApoFTO': 'ratio_log_ApoFTO',
            'ratio_TadFTO': 'ratio_log_TadFTO'
        }
    },
    'score3': {
        'primary': {
            'SKOV3_Apo_RIPK1_BT1_UT_zscore': 'zscore_RIPK1',
            'SKOV3_Apo_FTO_BT1_UT_zscore': 'zscore_ApoFTO',
            'SKOV3_Tad_FTO_BT1_UT_zscore': 'zscore_TadFTO'
        },
        'alternative': {
            'RIPK1_zscore': 'zscore_RIPK1',
            'ApoFTO_zscore': 'zscore_ApoFTO',
            'TadFTO_zscore': 'zscore_TadFTO'
        }
    },
    'score4': {
        'primary': {
            'SKOV3_Apo_RIPK1_BT1_UT_zscore': 'zscore_RIPK1',
            'SKOV3_Apo_FTO_BT1_UT_zscore': 'zscore_ApoFTO',
            'SKOV3_Tad_FTO_BT1_UT_zscore': 'zscore_TadFTO'
        },
        'alternative': {
            'RIPK1_zscore': 'zscore_RIPK1',
            'ApoFTO_zscore': 'zscore_ApoFTO',
            'TadFTO_zscore': 'zscore_TadFTO'
        }
    }
}

def get_column_mapping(df: pd.DataFrame, score: str, dataset: str) -> dict:
    """Determine which column mapping to use based on dataset and available columns."""
    if score not in score_column_mappings:
        return {}
        
    # Get mappings for this score
    mappings = score_column_mappings[score]
    
    primary = mappings['primary']
    if all(col in df.columns for col in primary):
        return primary
        
    alt = mappings['alternative']
    if all(col in df.columns for col in alt):
        return alt
        
    # Special handling for different datasets
    if dataset in ['mincov10_minedit0.05', 'mincov10_minedit0.05_sig']:
        # Try simpler column names
        if score in ['score1', 'score2']:
            simple_cols = {
                'ratio_RIPK1': 'ratio_RIPK1',
                'ratio_ApoFTO': 'ratio_ApoFTO',
                'ratio_TadFTO': 'ratio_TadFTO'
            }
            if all(col in df.columns for col in simple_cols):
                return simple_cols
        else:  # score3 and score4
            zscore_cols = {
                'RIPK1_zscore': 'zscore_RIPK1',
                'ApoFTO_zscore': 'zscore_ApoFTO',
                'TadFTO_zscore': 'zscore_TadFTO'
            }
            if all(col in df.columns for col in zscore_cols):
                return zscore_cols
    
    logging.warning(f"Dataset: {dataset}, Score: {score}")
    logging.warning(f"Available columns: {', '.join(df.columns)}")
    logging.warning(f"Expected columns not found")
    return {}

def load_common_sites(file_path: str = 'common_sites.csv') -> pd.DataFrame:
    """Load common sites with coordinates."""
    try:
        common_sites = pd.read_csv(file_path)
        # Create site identifier
        common_sites['site'] = common_sites.apply(
            lambda row: f"{row['chr']}_{row['start']}_{row['end']}", 
            axis=1
        )
        logging.info(f"Loaded {len(common_sites)} common sites")
        return common_sites
    except Exception as e:
        logging.error(f"Error loading common sites: {str(e)}")
        return None

def get_site_data(file_path: str, common_sites: pd.DataFrame) -> pd.DataFrame:
    """Extract data for common sites from CSV file."""
    try:
        # Read file without index for files like score_count_v1.csv
        df = pd.read_csv(file_path)
        
        if 'merged_column' not in df.columns:
            # Try reading with index for other format
            df = pd.read_csv(file_path, index_col=0)
            if 'merged_column' not in df.columns:
                logging.error(f"No merged_column found in {Path(file_path).name}")
                return None
        
        def extract_coords(x):
            try:
                parts = x.split(';')
                return f"{parts[0]}_{parts[1]}_{parts[2]}"
            except:
                return None
                
        df['temp_site'] = df['merged_column'].apply(extract_coords)
        
        site_ids = common_sites.apply(
            lambda row: f"{row['chr']}_{row['start']}_{row['end']}", 
            axis=1
        ).tolist()
        
        common_data = df[df['temp_site'].isin(site_ids)]
        
        if len(common_data) == 0:
            logging.error(f"No matching sites found in {Path(file_path).name}")
            return None
            
        if len(common_data) != len(common_sites):
            logging.warning(f"Found {len(common_data)}/{len(common_sites)} sites in {Path(file_path).name}")
        
        return common_data
        
    except Exception as e:
        logging.error(f"Error processing {Path(file_path).name}: {str(e)}")
        logging.error(f"At line: {e.__traceback__.tb_lineno}")
        return None

def plot_correlations(input_files: dict, common_sites: pd.DataFrame) -> None:
    """Generate correlation plots for each dataset and score type."""
    
    output_dir = Path('correlation_plots')
    output_dir.mkdir(exist_ok=True)
    
    for dataset, scores in input_files.items():
        logging.info(f"\nProcessing {dataset}")
        
        for score, metrics in scores.items():
            data_dict = {}
            
            for metric, file_path in metrics.items():
                common_data = get_site_data(file_path, common_sites)
                if common_data is not None:
                    # Get appropriate column mapping for this dataset
                    col_mapping = get_column_mapping(common_data, score, dataset)
                    
                    if col_mapping:
                        for orig_col, new_col in col_mapping.items():
                            if orig_col in common_data.columns:
                                data_dict[f"{new_col}_{metric}"] = common_data[orig_col]
            
            if data_dict:
                try:
                    # Create correlation matrix
                    corr_df = pd.DataFrame(data_dict).corr()
                    
                    plt.figure(figsize=(12, 10))
                    mask = np.triu(np.ones_like(corr_df), k=1)
                    
                    sns.heatmap(
                        corr_df,
                        mask=mask,
                        annot=True,
                        fmt='.2f',
                        cmap='coolwarm',
                        vmin=-1,
                        vmax=1,
                        center=0,
                        square=True,
                        annot_kws={'size': 8}
                    )
                    
                    plt.title(f'Correlation Matrix - {dataset} - {score}\n(n={len(common_sites)} sites)')
                    plt.tight_layout()
                    
                    output_base = output_dir / f'{dataset}_{score}'
                    plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    corr_df.to_csv(f'{output_base}_correlations.csv')
                    
                    logging.info(f"Saved correlation plot and data for {dataset} {score}")
                except Exception as e:
                    logging.error(f"Error creating plot for {dataset} {score}: {str(e)}")
            else:
                logging.error(f"No valid data found for {dataset} {score}")

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    common_sites = pd.read_csv('common_sites.csv')
    logging.info(f"Loaded {len(common_sites)} common sites")
    
    plot_correlations(input_files, common_sites)
