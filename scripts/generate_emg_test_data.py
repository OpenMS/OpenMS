#!/usr/bin/env python3
"""
Generate EMG (Exponentially Modified Gaussian) test data to reproduce issue #6239.

This script generates test data where:
1. A peak is positioned near the edge of the data range
2. The data is truncated to simulate a chromatogram ending abruptly
3. Without boundary constraints, the optimizer would fit RT outside the data range

The generated data can be used to verify the fix in EmgFitter1D.
"""

import math

def emg_function(t_values, height, width, symmetry, retention):
    """
    Simplified EMG function (doi=10.1.1.915.3568) Equation 9

    Parameters:
    - t_values: time points (list)
    - height: peak height (h)
    - width: peak width (w)
    - symmetry: peak symmetry/tau (s)
    - retention: retention time (z)
    """
    sqrt2pi = math.sqrt(2.0 * math.pi)
    emg_const = 2.4055
    sqrt_2 = math.sqrt(2.0)
    c = -emg_const / sqrt_2

    prefix = (height * width / symmetry) * sqrt2pi
    part1 = width * width / (2 * symmetry * symmetry)
    part2 = width / symmetry

    result = []
    for ti in t_values:
        diff = ti - retention
        exp_arg1 = part1 - (diff / symmetry)
        exp_arg2 = c * ((diff / width) - part2)

        # Avoid overflow
        if exp_arg1 > 700:
            yi = 0.0
        elif exp_arg2 > 700:
            yi = 0.0
        else:
            yi = prefix * math.exp(exp_arg1) / (1 + math.exp(exp_arg2))
        result.append(yi)

    return result


def generate_truncated_peak_data():
    """
    Generate test data for a peak near the upper edge that gets truncated.
    This reproduces the scenario from issue #6239.
    """
    # Peak parameters - peak center at RT=195, close to the truncation point at 200
    height = 50000.0
    width = 5.0
    symmetry = 5.0
    retention = 195.0

    # Full data range would be 100-300, but we truncate at 200
    t_full = [100.0 + i for i in range(200)]  # 100 to 299
    y_full = emg_function(t_full, height, width, symmetry, retention)

    # Truncate at RT=200
    truncation_point = 200.0
    t_truncated = [t for t in t_full if t <= truncation_point]
    y_truncated = [y_full[i] for i, t in enumerate(t_full) if t <= truncation_point]

    print("=" * 60)
    print("UPPER BOUND TEST DATA (Peak near upper edge)")
    print("=" * 60)
    print(f"True parameters: height={height}, width={width}, symmetry={symmetry}, retention={retention}")
    print(f"Data range: [{t_truncated[0]}, {t_truncated[-1]}]")
    print(f"Number of points: {len(t_truncated)}")
    print(f"Peak is at RT={retention}, data ends at RT={truncation_point}")
    print(f"Peak center is {truncation_point - retention} units from the upper edge")
    print()

    # Find the maximum intensity point in truncated data
    max_val = max(y_truncated)
    max_idx = y_truncated.index(max_val)
    print(f"Max intensity in truncated data: {max_val:.2f} at RT={t_truncated[max_idx]}")

    # Print C++ test code
    print("\n// C++ test data for upper bound test:")
    print("EmgModel::SamplesType truncated_samples;")
    for i in range(0, len(t_truncated), 10):  # Sample every 10th point for brevity
        print(f"truncated_samples.push_back(Peak1D({t_truncated[i]:.1f}, {y_truncated[i]:.2f}));")
    print("// ... (use interpolation step 1.0 to get all points)")

    return t_truncated, y_truncated


def generate_lower_bound_test_data():
    """
    Generate test data for a peak near the lower edge.
    """
    height = 50000.0
    width = 5.0
    symmetry = 5.0
    retention = 10.0  # Peak near lower edge

    # Full data range would be 0-200, but we truncate at lower end at 5
    t_full = [float(i) for i in range(200)]  # 0 to 199
    y_full = emg_function(t_full, height, width, symmetry, retention)

    # Truncate at RT=5 (lower end)
    truncation_point = 5.0
    t_truncated = [t for t in t_full if t >= truncation_point]
    y_truncated = [y_full[i] for i, t in enumerate(t_full) if t >= truncation_point]

    print("=" * 60)
    print("LOWER BOUND TEST DATA (Peak near lower edge)")
    print("=" * 60)
    print(f"True parameters: height={height}, width={width}, symmetry={symmetry}, retention={retention}")
    print(f"Data range: [{t_truncated[0]}, {t_truncated[-1]}]")
    print(f"Number of points: {len(t_truncated)}")
    print(f"Peak is at RT={retention}, data starts at RT={truncation_point}")
    print(f"Peak center is {retention - truncation_point} units from the lower edge")

    return t_truncated, y_truncated


def analyze_fitting_behavior():
    """
    Analyze what happens when fitting truncated data without bounds.

    The key insight: when data is truncated, the optimizer has no penalty
    for moving the peak center outside the data range, because there are
    no data points there to create residuals.
    """
    print("\n" + "=" * 60)
    print("ANALYSIS: Why RT can go outside range without bounds")
    print("=" * 60)

    # Parameters
    height = 50000.0
    width = 5.0
    symmetry = 5.0
    true_retention = 195.0

    # Data range
    t = [100.0 + i for i in range(101)]  # 100 to 200, truncated at 200
    y_data = emg_function(t, height, width, symmetry, true_retention)

    # Try fitting with different RT values
    test_retentions = [190.0, 195.0, 200.0, 205.0, 210.0, 220.0]

    print("\nResidual sum of squares for different RT guesses:")
    print("-" * 40)
    for test_rt in test_retentions:
        y_model = emg_function(t, height, width, symmetry, test_rt)
        rss = sum((ym - yd) ** 2 for ym, yd in zip(y_model, y_data))
        in_range = "IN RANGE" if 100 <= test_rt <= 200 else "OUTSIDE!"
        print(f"RT={test_rt:6.1f}: RSS = {rss:15.2f}  [{in_range}]")

    print()
    print("Note: Without boundary constraints, the optimizer may find that")
    print("RT values outside the data range give similar or even lower RSS,")
    print("because the model shape can shift without penalty from missing data points.")


def suggest_better_test_parameters():
    """
    Suggest test parameters that more reliably trigger the bug.
    """
    print("\n" + "=" * 60)
    print("SUGGESTED TEST PARAMETERS")
    print("=" * 60)

    print("""
For a more reliable test case, consider:

1. **More asymmetric peak** (larger symmetry parameter):
   - symmetry = 10.0 or higher
   - This creates a longer tail that extends beyond the truncation point

2. **Peak center closer to truncation**:
   - retention = 198.0 with truncation at 200.0
   - Only 2 units of buffer

3. **Wider peak**:
   - width = 10.0
   - More of the peak extends beyond truncation

Example parameters that should trigger the issue:
- height = 50000.0
- width = 8.0
- symmetry = 10.0
- retention = 198.0
- truncation = 200.0

With these parameters, a significant portion of the peak's tail is
cut off, which may cause the optimizer to shift RT beyond 200 to
better fit the visible portion of the peak.
""")


def analyze_extreme_case():
    """
    Analyze a more extreme case where the peak center is AT the truncation point.
    This is more likely to trigger the bug.
    """
    print("\n" + "=" * 60)
    print("EXTREME CASE: Peak center at truncation point")
    print("=" * 60)

    # Parameters - peak center exactly at truncation
    height = 50000.0
    width = 5.0
    symmetry = 5.0
    true_retention = 200.0  # Peak center AT truncation

    # Data truncated at 200
    t = [100.0 + i for i in range(101)]  # 100 to 200
    y_data = emg_function(t, height, width, symmetry, true_retention)

    print(f"True RT = {true_retention} (at truncation point)")
    print(f"Data range: [100, 200]")
    print()

    # Try fitting with different RT values
    test_retentions = [195.0, 198.0, 200.0, 202.0, 205.0, 210.0, 215.0, 220.0]

    print("Residual sum of squares for different RT guesses:")
    print("-" * 50)
    for test_rt in test_retentions:
        y_model = emg_function(t, height, width, symmetry, test_rt)
        rss = sum((ym - yd) ** 2 for ym, yd in zip(y_model, y_data))
        in_range = "IN RANGE" if 100 <= test_rt <= 200 else "OUTSIDE!"
        print(f"RT={test_rt:6.1f}: RSS = {rss:15.2f}  [{in_range}]")

    print()
    print("When the true peak is AT the truncation point, shifting RT")
    print("outside might not significantly increase the RSS because")
    print("half the peak is already missing from the data.")


def analyze_peak_beyond_truncation():
    """
    Analyze case where true peak center is BEYOND truncation point.
    This simulates a chromatogram that was cut off mid-peak.
    """
    print("\n" + "=" * 60)
    print("REAL BUG CASE: True peak beyond truncation")
    print("=" * 60)

    # Parameters - peak center BEYOND truncation
    height = 50000.0
    width = 5.0
    symmetry = 5.0
    true_retention = 210.0  # Peak center BEYOND truncation at 200

    # Data truncated at 200 - we only see the rising edge of the peak
    t = [100.0 + i for i in range(101)]  # 100 to 200
    y_data = emg_function(t, height, width, symmetry, true_retention)

    print(f"True RT = {true_retention} (BEYOND truncation at 200)")
    print(f"Data range: [100, 200]")
    print(f"We only see the rising edge of the peak!")
    print()

    # Find max in visible data
    max_val = max(y_data)
    max_idx = y_data.index(max_val)
    print(f"Max visible intensity: {max_val:.2f} at RT={t[max_idx]}")
    print()

    # What would an unconstrained optimizer find?
    test_retentions = [195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 230.0, 250.0]

    print("RSS for different RT guesses (optimizer might choose lowest):")
    print("-" * 50)
    min_rss = float('inf')
    best_rt = None
    for test_rt in test_retentions:
        y_model = emg_function(t, height, width, symmetry, test_rt)
        rss = sum((ym - yd) ** 2 for ym, yd in zip(y_model, y_data))
        in_range = "IN RANGE" if 100 <= test_rt <= 200 else "OUTSIDE!"
        marker = ""
        if rss < min_rss:
            min_rss = rss
            best_rt = test_rt
            marker = " <-- BEST SO FAR"
        print(f"RT={test_rt:6.1f}: RSS = {rss:15.2f}  [{in_range}]{marker}")

    print()
    print(f"An unconstrained optimizer would likely choose RT={best_rt}")
    if best_rt > 200:
        print("THIS IS OUTSIDE THE DATA RANGE - THE BUG!")
    print()
    print("The fix ensures the optimizer cannot return RT > 200 even")
    print("if the data suggests a peak beyond the truncation point.")


if __name__ == "__main__":
    generate_truncated_peak_data()
    print()
    generate_lower_bound_test_data()
    analyze_fitting_behavior()
    suggest_better_test_parameters()
    analyze_extreme_case()
    analyze_peak_beyond_truncation()
