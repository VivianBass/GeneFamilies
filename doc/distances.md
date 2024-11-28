

# ----------------------------------------------------------------------

# euclidean log2 Distances

A log₂ transformation applies the logarithm base 2 to your values. Here's what happens when you apply the transformation to values in the specified ranges:

---

### **1. Values Between 0 and 1**

- **Logarithmic Effect**: The logarithm of numbers between 0 and 1 is **negative** because the base (2) is greater than the number.
- **Transformation**: These values will transform into a range of negative numbers. For example:
  - log₂(0.5) = -1
  - log₂(0.25) = -2
  - log₂(0.1) ≈ -3.32
- **Behavior**: Values closer to 1 will transform to values closer to 0, while values closer to 0 will transform to large negative numbers.

---

### **2. Values Between 1 and 2**

- **Logarithmic Effect**: The logarithm of numbers between 1 and 2 is **positive but less than 1** because these numbers are greater than 1 but less than the base (2).
- **Transformation**: These values will transform into a range between 0 and 1. For example:
  - log₂(1) = 0
  - log₂(1.5) ≈ 0.585
  - log₂(2) = 1
- **Behavior**: The closer the value is to 2, the closer the transformed value will be to 1.

---

### **3. Values Between 2 and 10**

- **Logarithmic Effect**: The logarithm of numbers greater than 2 is **positive and greater than 1** because these values exceed the base.
- **Transformation**: These values will transform into a range of positive numbers greater than 1. For example:
  - log₂(2) = 1
  - log₂(4) = 2
  - log₂(8) = 3
  - log₂(10) ≈ 3.32
- **Behavior**: The larger the value, the greater the transformed value.

---

### **Summary of Effects**

| **Original Range** | **Transformed Range**          | **Key Features**                     |
| ------------------ | ------------------------------ | ------------------------------------ |
| 0 < x < 1          | Negative values (e.g., -1, -2) | Closer to 0 → closer to 0 post-log   |
| 1 ≤ x < 2          | 0 to 1                         | Closer to 2 → closer to 1 post-log   |
| 2 ≤ x ≤ 10         | 1 to ~3.32                     | Larger values → larger positive logs |

In short:

- Values between 0 and 1 shrink and become negative.
- Values between 1 and 2 transform to small positive numbers.
- Values greater than 2 scale to progressively larger positive values.


- add 1 so you dont have negative values !!!
- or Apply a pseudo-count (e.g., add a small constant like 1) to avoid taking the logarithm of zero.


# ---------------------------------------------------------------------


- cosine similarity, cosine distance, and angular distances
- Cosine similartity / Distance / Vectors, Scalars
- arccosine
- cosine Distance vs Cosine angle
- dot product of Vectors
- euclidean angles (complex numbers)

- result in angle (radian, euclidean angle (complex numbers ??)) or distance ???








---

### **Summary of Concepts and Explanation**

This explanation covers how to compute angular distances between vectors or scalar values using cosine similarity and arccosine. It compares the original Euclidean distance approach with angular measures to analyze directional relationships between gene expression data.

---

### **Key Concepts:**

#### **1. Cosine Similarity for Vectors**
- **Definition**: Measures how similar two vectors are by comparing their directions, independent of magnitude.

**Formula**:
$$
\text{cosine\ similarity}(u, v) = \frac{u \cdot v}{\|u\| \|v\|}
$$

Where:
- \( u \cdot v \): Dot product of vectors \( u \) and \( v \).
- \( \|u\|, \|v\| \): Vector magnitudes (norms).

**Angle Between Vectors**:
$$
\theta = \arccos(\text{cosine\_similarity}(u, v))
$$

- The result is an angle (\( \theta \)) in radians or degrees that quantifies directional similarity.

**Advantages**:
- Focuses on vector orientation, making it scale-invariant.
- Particularly useful for high-dimensional data (e.g., gene expression across tissues).

---

#### **2. Angular Distance Matrix**
- Compute pairwise cosine similarity for all vectors in a dataset.
- Convert similarities to angles using \( \arccos \).
- Store angles in a symmetric matrix, representing angular distances between all pairs.

---

#### **3. Comparison of Scalars**
- Scalars lack dimensionality, but you can adapt cosine-like methods to compare their direction (signs and magnitudes).

**Cosine-like Measure**: For two scalar values \( x1, x2 \):
$$
\text{Cosine\ Similarity}(x1, x2) = \frac{x1 \cdot x2}{|x1| \cdot |x2|}
$$

- Positive if \( x_1 \) and \( x_2 \) have the same sign.
- Negative if their signs differ.

**Angle from Scalars**:
$$
\text{angle} = \arccos(\text{Cosine\_Similarity})
$$

- This provides a measure of how aligned (or opposite) two scalar values are.




