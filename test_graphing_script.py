import matplotlib.pyplot as plt
import math

def log_migration(length, min_bp, max_bp, gel_height=10):
    """Map fragment length to gel position (log scale)."""
    log_min = math.log10(min_bp)
    log_max = math.log10(max_bp)
    norm_log = (math.log10(length) - log_min) / (log_max - log_min)
    return gel_height * (1 - norm_log)

def draw_gel_filtered(
    ladder,
    fragments,
    overlapping_fragments,
    min_len,
    max_len,
    resolution=0.05,
    gel_height=10
):
    # Filter based on cutoffs
    ladder = [f for f in ladder if min_len <= f <= max_len]
    fragments = [f for f in fragments if min_len <= f <= max_len]
    overlapping_fragments = [f for f in overlapping_fragments if min_len <= f <= max_len]

    fig, ax = plt.subplots(figsize=(2, 8))  # One ladder and one sample lane

    ax.set_xlim(0, 2)
    ax.set_ylim(0, gel_height)
    ax.set_facecolor("white")
    ax.set_ylabel("Fragment length (bp)", color='black')
    ax.tick_params(axis='y', colors='black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks([])

    # Invert y-axis (smaller bp at bottom)
    ax.invert_yaxis()

    # Y-axis ticks (log-spaced)
    tick_lengths = [100, 200, 500, 1000, 2000, 5000, 10000]
    tick_lengths = [t for t in tick_lengths if min_len <= t <= max_len]
    tick_positions = [log_migration(t, min_len, max_len, gel_height) for t in tick_lengths]
    ax.set_yticks(tick_positions)
    ax.set_yticklabels([str(t) for t in tick_lengths])

    # Draw ladder (lane 1)
    x_ladder = 0.5
    for frag in ladder:
        spread = resolution * frag
        y1 = log_migration(frag + spread, min_len, max_len, gel_height)
        y2 = log_migration(frag - spread, min_len, max_len, gel_height)
        ax.add_patch(plt.Rectangle((x_ladder - 0.35, y1), 0.7, y2 - y1, color='black'))

    ax.text(x_ladder, gel_height + 0.3, "Ladder", ha='center', color='black', fontsize=9)

    # Draw sample (lane 2)
    x_sample = 1.5
    for frag in fragments:
        spread = resolution * frag
        y1 = log_migration(frag + spread, min_len, max_len, gel_height)
        y2 = log_migration(frag - spread, min_len, max_len, gel_height)
        color = 'red' if frag in overlapping_fragments else 'gray'
        ax.add_patch(plt.Rectangle((x_sample - 0.35, y1), 0.7, y2 - y1, color=color))

    ax.text(x_sample, gel_height + 0.3, "Sample", ha='center', color='black', fontsize=9)

    plt.tight_layout()
    plt.show()

# ✅ Example call:
ladder = [100, 200, 400, 800, 1600, 3200, 6400]
fragments = [120, 250, 800, 900, 3100, 5000]
overlapping = [800, 3100]

draw_gel_filtered(
    ladder=ladder,
    fragments=fragments,
    overlapping_fragments=overlapping,
    min_len=100,
    max_len=5000,
    resolution=0.05  # 5% spread
)