# How the Van Allen Belt Flux is Calculated for Each Grid Point
# Used in generate_sample_data.py

This is a step by step explanation of the function
`van_allen_belt_model(x, y, z)`, describing how it computes the flux value for a single point in space.

---

## 1. Constants

- Earth radius: **Rₑ = 6371 km**
- Magnetic dipole tilt: **θ = 11°**
- Noise and background (applied at the end):
  - Multiply flux by \( 1 + 0.1\,\mathcal{N}(0,1) \) (clipped to ≥ 0.1)
  - Add a random offset \( \mathcal{U}(0, 10^4) \)

---

## 2. Rotate into Magnetic Coordinates

Incoming coordinates are geographic (rotation-axis aligned). They are rotated into magnetic coordinates:

\[
\begin{aligned}
x_{mag} &= x \cos\theta + z \sin\theta \\
y_{mag} &= y \\
z_{mag} &= -x \sin\theta + z \cos\theta
\end{aligned}
\]

---

## 3. Magnetic Geometry and L-Shell

From the rotated coordinates:

\[
\begin{aligned}
r_{mag} &= \sqrt{x_{mag}^2 + y_{mag}^2 + z_{mag}^2} \\
\lambda &= \arcsin\left(\frac{|z_{mag}|}{r_{mag}}\right) \\
L &= \frac{r_{mag}}{R_E \cos^2 \lambda}
\end{aligned}
\]

If \( r_{mag} \leq R_E \), the point is inside Earth → **flux = 0**.

---

## 4. Flux Contributions

Flux is computed as the sum of several Gaussian components:

### Inner Belt (L = 1.2–2.8)
\[
flux \;{+}{=} 2\times 10^7\; e^{-(L-1.6)^2 / (2(0.4)^2)} \; (1 + 2\sin^2 \lambda)
\]

### Slot Region (L = 2.5–3.2)
\[
flux \;{+}{=} 10^5 \; e^{-(L-2.8)^2 / (2(0.2)^2)}
\]

### Outer Belt (L = 3.0–8.0)
\[
flux \;{+}{=} 10^7 \; e^{-(L-4.5)^2 / (2(1.2)^2)} (1 + 1.5 \sin^2 \lambda)(1 + 0.3 \cos \phi)
\]

where \( \phi = \arctan2(y_{mag}, x_{mag}) \).

### Plasma Sheet Tail (L > 6)
\[
flux \;{+}{=} 5\times 10^5 \; e^{-(L-8)^2 / (2(2)^2)}
\]

---

## 5. Noise and Background

Finally, if flux > 0:

\[
\begin{aligned}
flux &\leftarrow flux \times \max(0.1, 1 + 0.1\,\mathcal{N}(0,1)) \\
flux &\leftarrow flux + \mathcal{U}(0, 10^4)
\end{aligned}
\]

This ensures every point has at least a small nonzero value.

---

## 6. Worked Example

Point: \( (x,y,z) = (2R_E, 0, 0) = (12742, 0, 0) \; km \)

**Step 1. Rotate into magnetic frame**

\[
\begin{aligned}
x_{mag} &= 12742 \cos 11^\circ \approx 12508 \\
y_{mag} &= 0 \\
z_{mag} &= -12742 \sin 11^\circ \approx -2431
\end{aligned}
\]

**Step 2. Geometry**

\[
r_{mag} = 12742, \quad
\lambda = \arcsin(2431/12742) = 11^\circ, \quad
L = 2.08
\]

**Step 3. Flux contributions**

- Inner belt is active:  
\( flux \approx 1.06 \times 10^7 \)  
- Slot region: none (L < 2.5)  
- Outer belt: none (L < 3)  
- Plasma sheet: none (L < 6)

**Step 4. Noise + background**

Apply random multiplier (≈1.003) and random offset (≈4389).  

Final flux ≈ **1.09 × 10⁷** (dimensionless intensity).

---

## 7. Important Note on Units

- The script uses arbitrary amplitudes (so no reference spectra, no
calibrations, and arbitrary units).  
- Output is for visualization purposes only.

