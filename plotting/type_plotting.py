import matplotlib.pyplot as plt


def plot_parallel_performance(parallel_types, execution_times):
    """
    Create a bar chart comparing execution times for different parallel types

    Parameters:
    - parallel_types: List of parallel types (e.g., ['Serial', 'OMP', 'Thread', 'MPI'])
    - execution_times: List of corresponding execution times
    """
    plt.figure(figsize=(10, 6))

    # Create bar chart with multiple colors
    bars = plt.bar(parallel_types, execution_times,
                   color=['blue', 'green', 'red', 'purple'])

    plt.title('ACO Parallel Performance Comparison', fontsize=16)
    plt.xlabel('Parallel Type', fontsize=12)
    plt.ylabel('Execution Time (milliseconds)', fontsize=12)

    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:,.0f} ms',
                 ha='center', va='bottom', fontsize=10)

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save figure
    plt.savefig('../res/parallel_performance.png',
                dpi=300, bbox_inches='tight')
    print("Figure saved to ./results/parallel_performance.png")

    plt.show()


# Example usage
parallel_types = ['Serial', 'OMP', 'Thread_32', 'MPI']
execution_times = [41337, 3123, 2287, 2099]
plot_parallel_performance(parallel_types, execution_times)
