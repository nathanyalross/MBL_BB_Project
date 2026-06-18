import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from dataclasses import dataclass

def data_analysis(data_1, data_2):
        """
        Collection of data tests to compare two datasets. First, test for normalization then compare depending on outcome
        
        :param data_1: Dataset as a list of values
        :param data_2: Dataset as a list of values

        Returns: test statistic, p-value for comparison, and statistical test used
        """
        
        #First, test both datasets for normal distribution
        norm_data_1 = stats.shapiro(data_1)
        norm_data_2 = stats.shapiro(data_2)

        #Depending on outcomes, use either paired T-test (norm) or Wilcoxon signed-rank test (non-norm)
        if (norm_data_1.pvalue < 0.05) or (norm_data_2.pvalue < 0.05): #Select the p-value, a statistic is also returned here
                w_stat, p_value = stats.wilcoxon(data_1,data_2)
                return (w_stat, p_value, 'Wilcoxon signed-rank')
        else:
                t_stat, p_value = stats.ttest_rel(data_1,data_2)
                return (t_stat, p_value, 'paired T-test')

def plot_line(dataset,
                        color,
                        title, 
                        ylabel, 
                        xlabel='Sweeps',
                        drug_region=(5, 10),
                        drug_alpha=0.2,
                        drug_color='gray',
                        figsize=(8, 6),
                        show_legend=True,
                        normalization = False):
    """
    Plot a single dataset as a line with SEM (standard error of mean) shading.
    
    Parameters:
    -----------
    dataset: Cell data as a DataFrame with rows as sweeps and columns as cells
    color: hex color of plotted line
    title: Plot Title
    xlabel: x-axis label
    ylabel: y-axis label
    drug_region: tuple (start_idx, end_idx) for highlighting drug application region
    drug_alpha: transparency of drug region shading
    drug_color: color of drug region shading
    figsize: tuple (width, height) for figure size
    show_legend: whether to show legend
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    
    fig, ax = plt.subplots(figsize=figsize)

    #Optionally Normalize the data for each cell
    if normalization is not False:
        #Calculate baseline frequency by averaging first 5 sweeps
        bl_mean = dataset.iloc[0:4,:].mean(axis = 0)
        #Divide all values by bl mean
        data = dataset.div(bl_mean)
    else:
        data = dataset
    
    # Calculate mean and SEM GFP+
    mean = data.mean(axis=1)  # Mean across columns (cells)
    sem = data.std(axis=1) / np.sqrt(data.shape[1])
        
    x = np.arange(1, (len(mean)+1))
        
    # Plot mean line
    ax.plot(x, mean, color=color, linewidth=2)
        
    # Plot SEM shading
    ax.fill_between(x, mean - sem, mean + sem, alpha=0.3, color=color)

    # Highlight drug application region
    if drug_region is not None:
        ax.axvspan(drug_region[0], drug_region[1], 
                   alpha=drug_alpha, color=drug_color, 
                   label='Drug Application', zorder=0)
    
    # Customize axes
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add title if provided
    if title:
        ax.set_title(title, fontsize=14)
    
    # Add legend
    if show_legend:
        ax.legend(frameon=False)
    
    plt.tight_layout()
    
    return fig

def plot_bar_two_datasets(data_1, bar1_color, data_2, bar2_color, x_tick_labels, plot_title, y_axis_label):
        """
        Function to make a simple bar graph to compare two datasets
        
        :param data_1: List of datapoints for dataset 1
        :param bar1_color: Hexcode color for first bar 
        :param data_2: List of datapoints for dataset 2
        :param bar2_color: Hexcode color for second bar 
        :param x_tick_labels: List of two x-tick labels, must be a list of strings
        :param plot_title: Title of plot as a string
        :param y_axis_label: y-axis label for graph

        Returns: matplotlib figure object
        """
        # Perform data analysis
        test_stat, p_value, test_type = data_analysis(data_1, data_2)


        # Bar edge properties
        bar_edge_color = 'black'
        bar_edge_width = 1

        # Font sizes
        title_fontsize = 18
        axis_label_fontsize = 16

        # Bar width
        bar_width = 0.6

        # Figure size (width, height in inches)
        fig_width = 4
        fig_height = 6

        # ==============================================================================
        # CREATE THE PLOT
        # ==============================================================================

        # Calculate means for bar heights
        data_1_mean = np.mean(data_1)
        data_2_mean = np.mean(data_2)

        data_1_sem = np.std(data_1) / np.sqrt(len(data_1))
        data_2_sem = np.std(data_2) / np.sqrt(len(data_2))

        # Set up the figure with high DPI for publication quality
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

        # X-axis positions
        x_positions = np.arange(2)  # Two bars

        # Create bars
        bars = ax.bar(x_positions, [data_1_mean, data_2_mean], bar_width,
                        color=[bar1_color, bar2_color],
                        edgecolor=bar_edge_color,
                        linewidth=bar_edge_width,
                        yerr=[data_1_sem, data_2_sem],
                        capsize=5,
                        error_kw={'linewidth': 2})

        # Add individual data points with jitter
        np.random.seed(42)  # For reproducible jitter

        # Plot GFP+ data points
        x_jitter_gfp_pos = np.random.normal(0, 0.05, size=len(data_1))
        ax.scatter(x_positions[0] + x_jitter_gfp_pos, data_1,
                color='#FFFFFF', s=50, alpha=0.6,
                edgecolors='black', linewidths=0.5, zorder=3)

        # Plot GFP- data points
        x_jitter_gfp_neg = np.random.normal(0, 0.05, size=len(data_2))
        ax.scatter(x_positions[1] + x_jitter_gfp_neg, data_2,
                color='#FFFFFF', s=50, alpha=0.6,
                edgecolors='black', linewidths=0.5, zorder=3)
        # Add paired connecting lines
        for i in range(min(len(data_1), len(data_2))):
            ax.plot([x_positions[0], x_positions[1]],
                    [data_1[i], data_2[i]],
                    color='gray', alpha=0.5, linewidth=1, zorder=2)

        # Set labels and title
        ax.set_ylabel(y_axis_label, fontsize=axis_label_fontsize, fontweight='bold')
        ax.set_title(plot_title, fontsize=title_fontsize, fontweight='bold', pad=20)

        # Set x-axis ticks and labels
        ax.set_xticks(x_positions)
        ax.set_xticklabels(x_tick_labels, fontsize=16, fontweight = 'bold')

        # Set y-axis tick label size
        ax.tick_params(axis='y', labelsize=12)

        # Remove top and right spines (borders)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Make left and bottom spines thicker
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)

        # Add a line at 0
        ax.axhline(y=0, color='black', linewidth=1)

        # Determine appropriate y-axis limits
        all_data = np.concatenate([data_1, data_2])
        y_min = min(0, np.min(all_data), data_1_mean - data_1_sem, data_2_mean - data_2_sem)
        y_max = max(0, np.max(all_data), data_1_mean + data_1_sem, data_2_mean + data_2_sem)

        # Add some padding (10% of the range)
        y_range = y_max - y_min
        y_min -= 0.1 * y_range
        y_max += 0.1 * y_range

        ax.set_ylim(y_min, y_max)

        #Add a marker for significance if p-value is less than 0.05
        if p_value < 0.05:
                # Determine significance level
                if p_value < 0.0001:
                        sig_marker = '****'
                elif p_value < 0.001:
                        sig_marker = '***'
                elif p_value < 0.01:
                        sig_marker = '**'
                else:
                        sig_marker = '*'
                
                # Position the line near the top of the plot
                # Calculate this before setting y_lim
                all_data = np.concatenate([data_1, data_2])
                y_min_temp = min(0, np.min(all_data), data_1_mean - data_1_sem, data_2_mean - data_2_sem)
                y_max_temp = max(0, np.max(all_data), data_1_mean + data_1_sem, data_2_mean + data_2_sem)
                y_range_temp = y_max_temp - y_min_temp
                
                # Position line at ~90% of the y-axis height
                bar_height = y_max_temp + 0.08 * y_range_temp
                
                # Draw just the horizontal significance line
                ax.plot([x_positions[0], x_positions[1]], [bar_height, bar_height], 
                        'k-', linewidth=1.5)
                
                # Add the asterisk(s) above the line
                ax.text((x_positions[0] + x_positions[1]) / 2, bar_height + 0.01 * y_range_temp, 
                        sig_marker, ha='center', va='bottom', fontsize=16, fontweight='bold')

        # Adjust layout to prevent label cutoff
        fig.tight_layout()

        # Return the figure object
        return fig, test_stat, p_value, test_type

@dataclass
class CumulativeProbabilityStats:
    """Statistical results from the cumulative probability comparison."""
    ks_statistic: float
    ks_pvalue: float
    significance_label: str
    baseline_n: int
    drug_n: int
    baseline_median: float
    drug_median: float
    baseline_mean: float
    drug_mean: float

    def __repr__(self):
        return (
            f"KS Statistic:       {self.ks_statistic:.4f}\n"
            f"KS p-value:         {self.ks_pvalue:.4f}\n"
            f"Significance:       {self.significance_label}\n"
            f"Baseline  — n={self.baseline_n}, mean={self.baseline_mean:.3f}, median={self.baseline_median:.3f}\n"
            f"Drug      — n={self.drug_n},  mean={self.drug_mean:.3f}, median={self.drug_median:.3f}\n"
        )

def _significance_label(p: float) -> str:
    """Convert p-value to asterisk notation."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

def plot_cumulative_probability(
    df: pd.DataFrame | list[pd.DataFrame],  # accept either
    column: str,
    n_bootstrap: int = 1000,
    confidence: float = 0.95,
    baseline_sweeps: tuple = (0, 5),
    drug_sweeps: tuple = (5, 10),
    ax: plt.Axes = None,
    title: str = None,
) -> tuple[plt.Figure, plt.Axes, CumulativeProbabilityStats]:
    """
    Plot cumulative probability distributions with bootstrapped confidence intervals
    and a two-sample KS test for baseline vs. drug application periods.

    Parameters
    ----------
    df              : DataFrame containing 'Sweep Number' and the target column
    column          : Column name to plot (e.g. 'Peak (pA)', 'AUC', 'Rise Time (ms)')
    n_bootstrap     : Number of bootstrap iterations (default 1000)
    confidence      : Confidence interval width (default 0.95 → 95% CI)
    baseline_sweeps : Inclusive sweep range for baseline (default (1, 5))
    drug_sweeps     : Inclusive sweep range for drug (default (7, 11))
    ax              : Optional existing Axes to plot onto
    title           : Optional plot title override

    Returns
    -------
    fig             : Matplotlib Figure
    ax              : Matplotlib Axes
    stats           : CumulativeProbabilityStats dataclass with all statistical results
    """

    # ── 1. Normalize input and extract data ────────────────────────────────────────────
    if isinstance(df, pd.DataFrame):
        df_list = [df]
    else:
        df_list = df

    # Then aggregate all data across DataFrames before proceeding
    import numpy as np
    baseline_data = np.concatenate([
        d[d["Sweep Number"].between(*baseline_sweeps)][column].dropna().values
        for d in df_list
    ])
    drug_data = np.concatenate([
        d[d["Sweep Number"].between(*drug_sweeps)][column].dropna().values
        for d in df_list
    ])

    if len(baseline_data) == 0 or len(drug_data) == 0:
        raise ValueError(
            f"No data found for column '{column}' in one or both sweep ranges. "
            f"Check your sweep numbers and column name."
        )

    # ── 2. KS Test ────────────────────────────────────────────────────────────
    ks_stat, ks_p = stats.ks_2samp(baseline_data, drug_data)
    sig_label = _significance_label(ks_p)

    result_stats = CumulativeProbabilityStats(
        ks_statistic    = ks_stat,
        ks_pvalue       = ks_p,
        significance_label = sig_label,
        baseline_n      = len(baseline_data),
        drug_n          = len(drug_data),
        baseline_median = float(np.median(baseline_data)),
        drug_median     = float(np.median(drug_data)),
        baseline_mean   = float(np.mean(baseline_data)),
        drug_mean       = float(np.mean(drug_data)),
    )

    # ── 3. Bootstrap resampling ───────────────────────────────────────────────
    alpha = 1 - confidence

    def bootstrap_cdf(data, n_boot):
        x = np.sort(data)
        boot_cdfs = np.zeros((n_boot, len(x)))
        for i in range(n_boot):
            resample = np.random.choice(data, size=len(data), replace=True)
            boot_cdfs[i] = np.searchsorted(np.sort(resample), x, side="right") / len(resample)
        median_cdf = np.median(boot_cdfs, axis=0)
        lower_ci   = np.percentile(boot_cdfs, 100 * (alpha / 2), axis=0)
        upper_ci   = np.percentile(boot_cdfs, 100 * (1 - alpha / 2), axis=0)
        return x, median_cdf, lower_ci, upper_ci

    b_x, b_med, b_lo, b_hi = bootstrap_cdf(baseline_data, n_bootstrap)
    d_x, d_med, d_lo, d_hi = bootstrap_cdf(drug_data,     n_bootstrap)

    # ── 4. Plot ───────────────────────────────────────────────────────────────
    BASELINE_COLOR = "#f0c6c6"
    DRUG_COLOR     = "#a12c2a"

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.get_figure()

    ax.fill_between(b_x, b_lo, b_hi, alpha=0.20, color=BASELINE_COLOR, linewidth=0)
    ax.fill_between(d_x, d_lo, d_hi, alpha=0.20, color=DRUG_COLOR,     linewidth=0)
    ax.step(b_x, b_med, where="post", color=BASELINE_COLOR, linewidth=2, label="Baseline")
    ax.step(d_x, d_med, where="post", color=DRUG_COLOR,     linewidth=2, label="Drug")

    # ── 5. KS significance annotation ────────────────────────────────────────
    # Place bracket + asterisks at the top of the plot
    ax_top = 0.97
    bracket_y = ax_top - 0.04

    # Find x positions of each curve's median (for bracket anchoring)
    x_left  = float(np.percentile(baseline_data, 50))
    x_right = float(np.percentile(drug_data, 50))
    x_mid   = (x_left + x_right) / 2

    # Draw bracket in axes-fraction coordinates for robustness
    ax.annotate(
        "", xy=(x_right, bracket_y), xytext=(x_left, bracket_y),
        xycoords=("data", "axes fraction"),
        textcoords=("data", "axes fraction"),
        arrowprops=dict(arrowstyle="-", color="black", lw=1.2,
                        connectionstyle="bar,fraction=0.15"),
    )
    ax.text(
        x_mid, ax_top,
        sig_label,
        ha="center", va="bottom",
        fontsize=14 if sig_label != "ns" else 11,
        fontweight="bold" if sig_label != "ns" else "normal",
        transform=ax.get_xaxis_transform(),  # x=data, y=axes fraction
    )

    # ── 6. Labels & formatting ────────────────────────────────────────────────
    ax.set_xlabel(column, fontsize=12)
    ax.set_ylabel("Cumulative Probability", fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_title(
        title or f"Cumulative Probability — {column}",
        fontsize=13, fontweight="bold"
    )

    p_display = f"p < 0.001" if ks_p < 0.001 else f"p = {ks_p:.3f}"
    ci_patch = Patch(color="gray", alpha=0.3,
                     label=f"{int(confidence*100)}% CI (bootstrap, n={n_bootstrap})")
    ks_patch = Patch(color="none",
                     label=f"KS = {ks_stat:.3f}, {p_display}")

    ax.legend(handles=[ax.lines[0], ax.lines[1], ci_patch, ks_patch], fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)
    fig.tight_layout()

    return fig, ax, result_stats