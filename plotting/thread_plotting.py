import matplotlib.pyplot as plt
import os


def plot_thread_performance(threads, times, output_dir='./'):
    """
    Create a line plot of performance vs number of threads with equidistant x-axis

    Parameters:
    - threads: List of thread counts
    - times: List of corresponding execution times
    - output_dir: Directory to save the figure
    """
    plt.figure(figsize=(12, 7))

    # Plot the line graph with equidistant x-axis points
    plt.plot(range(len(threads)), times, marker='o',
             linestyle='-', linewidth=2, markersize=8)

    # Set x-ticks to match thread numbers
    plt.xticks(range(len(threads)), [str(t) for t in threads])

    plt.title('GA a280.tsp Performance: MPI vs Execution Time', fontsize=16)
    plt.xlabel('Number of Threads', fontsize=12)
    plt.ylabel('Execution Time (milliseconds)', fontsize=12)

    plt.grid(True, linestyle='--', alpha=0.7)

    # Annotate each point with its time
    for i, (thread, time) in enumerate(zip(threads, times)):
        plt.annotate(f'{time} ms',
                     (i, time),
                     xytext=(0, 10),
                     textcoords='offset points',
                     ha='center',
                     fontsize=9)

    # Save figure
    output_path = os.path.join(output_dir, 'thread_performance.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to {output_path}")

    plt.tight_layout()
    plt.show()


# Example usage
threads = [2, 4, 8, 16, 32, 64]
times = [161058, 81811, 42996, 23705, 14370, 36908]
plot_thread_performance(threads, times, output_dir='../res')
