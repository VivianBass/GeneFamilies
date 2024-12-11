
## **Statistical Analysis Methods**

### 1. T-Test Analysis

- **Purpose**: Compare means between groups (orthologs vs. paralogs)
- **Requirements**:
  - Normal distribution
  - Equal variances
- **Hypotheses**:
  - H₀: Group means are equal
  - H₁: Group means differ

### 2. Wilcoxon Rank-Sum Test

- **Purpose**: Non-parametric alternative to t-test
- **Use Case**: When normality assumptions aren't met
- **Approach**: Compares rank distributions
- **Effect Size**: 
  - Small: ≈0.1
  - Medium: ≈0.3
  - Large: ≈0.5

### 3. Implementation Notes

- **Test Direction**:
  - Use "two.sided" for unbiased comparison
  - Directional tests ("greater"/"less") affect interpretation
- **Group Order**:
  - Matters for one-sided tests
  - Irrelevant for two-sided tests

### 4. Multiple Testing Correction

- **Method**: Benjamini-Hochberg (BH)
- **Purpose**: Control False Discovery Rate (FDR)
- **Implementation**:
  ```R
  p.adjust(p_values, method = "BH")
  ```

### 5. Key Considerations

- Choose test based on data distribution
- Consider sample size and variance
- Apply appropriate multiple testing correction
- Interpret results in biological context






# ---------------------------------------------------------------------


### Detailed Summary: T-Test Analysis and P-Value Adjustments

#### **1. T-Test Basics**
- **Purpose**: Compare the means of two groups to determine if they differ significantly.
- **Hypotheses**:
  - \( H_0 \): The means of the two groups are equal (\( \mu_1 = \mu_2 \)).
  - \( H_1 \): One mean is greater than the other (\( \mu_1 > \mu_2 \), for one-sided tests).
- **T-Statistic Formula**:
  \[
  t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}
  \]
  - \( \bar{x}_1, \bar{x}_2 \): Sample means for the two groups.
  - \( s_1^2, s_2^2 \): Variances of the two groups.
  - \( n_1, n_2 \): Sample sizes for the two groups.
  - Measures the difference in means relative to variability. Larger \( t \)-values indicate greater differences.

#### **2. Degrees of Freedom (df)**
- Determines the shape of the t-distribution used for p-value calculations.
- For unpaired t-tests:
  \[
  df = n_1 + n_2 - 2
  \]

#### **3. P-Value Calculation**
- P-values represent the probability of observing a \( t \)-statistic as extreme as the observed one, assuming \( H_0 \) is true.
- **One-Sided Test (alternative = "greater")**:
  \[
  p = P(T \geq t_{\text{observed}})
  \]
  Where \( T \) follows a t-distribution with \( df \) degrees of freedom.
- For a two-sided test, \( p \) includes both tails of the distribution.

#### **4. Multiple Comparisons and Adjusted P-Values**
When performing multiple tests, the likelihood of false positives increases. To address this:
- **Benjamini-Hochberg (BH) Adjustment** controls the False Discovery Rate (FDR):
  1. Sort raw p-values (\( p_1, p_2, \dots, p_m \)) in ascending order.
  2. Compute adjusted p-values:
     \[
     p_{\text{adjusted}, i} = \min\left(\frac{p_i \cdot m}{i}, 1\right)
     \]
     Where \( m \) is the total number of tests, and \( i \) is the rank of the p-value.
  3. Ensures fewer false positives while maintaining statistical power.

#### **5. Interpreting Results**
- **Significance Thresholds**:
  - \( p < 0.001 \): Highly significant (***)
  - \( p < 0.01 \): Very significant (**)
  - \( p < 0.05 \): Significant (*)
  - \( p \geq 0.05 \): Not significant (ns).
- **Negative t-Values**:
  - Indicate that the second group’s mean is larger than the first. Sign flips if the group order changes, but p-values remain unchanged.
  
#### **6. rstatix T-Test Output**
The `t_test` function provides detailed output for each pairwise comparison:
- `.y.`: Dependent variable being tested (e.g., `Distance`).
- `group1` and `group2`: The two groups being compared.
- `n1` and `n2`: Sample sizes of the groups.
- `statistic`: T-statistic value.
- `df`: Degrees of freedom.
- `p`: Raw p-value.
- `p.adj`: Adjusted p-value using BH correction.
- `p.adj.signif`: Significance annotation (***, **, *, ns).
- `test_type`: Specifies the test performed (e.g., t-test).

This detailed output provides the statistical evidence needed to assess differences while controlling for multiple testing.