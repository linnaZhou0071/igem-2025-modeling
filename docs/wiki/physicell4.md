# PhysiCell4 — wiki manuscript

> Converted from the locally retained source document `PhysiCell4.docx`. Edited for portfolio readability: internal layout notes and conversion artifacts have been removed. The [published team wiki](https://2025.igem.wiki/sjtu-biox-shanghai/model) presents the final results; see [editorial notes](../EDITORIAL_NOTES.md) before interpreting its claims.

### Results-Validation of the Planting and Treatment Capacity

**2D-Ocurrence, Recurrence and Treatment Models**

- First, by comparing the ***Occurrence* and *Recurrence* models**, we validated the engineered strain’s **therapeutic efficacy.**

[Fig.26: Occurrence model simulation](https://static.igem.wiki/teams/5615/wiki/model/2d-ad-trial2-1.webp)

[Fig.27: Recurrence model simulation](https://static.igem.wiki/teams/5615/wiki/model/2d-recurrence-1.webp)

In the ***Occurrence* model**, when no engineered *S. epidermidis* is present, ***S. aureus* microcolonies expand** continuously and show signs of **biofilm-mediated rupture and invasion**. **T-cell cytokine** concentrations also remain **high**, indicating an inflamed microenvironment.

In contrast, in the **recurrence model**, the **engineered *S. epidermidis* persistently suppresses *S. aureus* growth**, resulting in both populations remaining near steady state. Consequently, **T-cell cytokine levels stay low**, indicating effective disease control.

![Fig.28 Changes in population plot and substance concentration in the Occurrence Model](media/physicell4/media/image1.gif)

Fig.28 Changes in population plot and substance concentration in the Occurrence Model

![Fig.29 Changes in population plot and substance concentration in the Recurrence Model](media/physicell4/media/image2.gif)

Fig.29 Changes in population plot and substance concentration in the Recurrence Model

As shown, population plots of cell counts again confirm suppression by the engineered strain. **Field simulations** show that in ***Occurrence***, **toxin concentration rises** with increasing *S. aureus* abundance (visible after ~40 h); in ***Recurrence***, toxin remains below visible levels, **cytokine levels decline** overall, and show clear **spatial correlation with VHH**.

In the *Occurrence* spatial fields we also observe successful **biofilm modeling**: after EPS-defined biofilms form, AMP is rapidly trapped and the spatial extent of bacteriocin_bound expands.

- Next, we used the ***Treatment* model** to test the engineered strain’s ability to **occupy the niche** and establish **long-term colonization**.

[Fig.30: Treatment model simulation](https://static.igem.wiki/teams/5615/wiki/model/2d-treatment-1.webp)

The simulation movie shows that during an **AD flare** when ***S. aureus* dominates**, environmental resources are heavily consumed by *S. aureus* and **cytokine levels are high**. **Massive seeding of *S. epidermidis* instantly reduces cytokine** concentrations, but may **exceed** the system’s **carrying capacity**, causing large-scale microbial mortality and a transition to a new steady state. After stabilization, the **engineered strain persists**, although cytokine levels show some oscillation.

![Fig.31 Changes in population plot and substance concentration in the Treatment Model](media/physicell4/media/image3.gif)

Fig.31 Changes in population plot and substance concentration in the Treatment Model

As shown, after colonization **engineered *S. epidermidis* secretes abundant VHH**, rapidly lowering T-cell cytokine and toxin concentrations and achieving successful establishment. However, **around ~60 h post-colonization cytokine levels rebound while VHH declines**, indicating that colonization alone may not sustain therapeutic control in severe AD. Therefore for **severe AD** we still recommend **subcutaneous antibody injection** combined with our **engineered strain as adjuvant therapy to prevent relapse**.

**3D-Ocurrence and Recurrence Models**

As the 2D results suggested that the engineered strain was **better suited to preventing relapse than treating severe AD**, we **only** used 3D simulations to study the strain’s role in **recurrence**.

Fig.32(a)/(b) 3D Simulation Movie with Th2 cytokine; left, Occurrence 3D; right, Recurrence 3D. [Paired simulation](https://static.igem.wiki/teams/5615/wiki/model/paired3d-ad-cytokine.webp) · [Companion panel](https://static.igem.wiki/teams/5615/wiki/model/33b-new.webp)

Fig.33(a)/(b) 3D Simulation Movie with toxin; left, Occurrence 3D; right, Recurrence 3D. [Paired simulation](https://static.igem.wiki/teams/5615/wiki/model/paired3d-ad-toxin.webp) · [Companion panel](https://static.igem.wiki/teams/5615/wiki/model/34b-new.webp)

In terms of concentration fields, **cytokine and toxin ranges** in the ***Recurrence*** 3D model (with engineered strain) are **lower** than in the *Occurrence* 3D model. Spatially, **T-cell cytokines** tend to **cluster in the dermis**, whereas **toxin** shows **localized bursts within *S. aureus*** populations at the **skin surface**; *S. aureus* colonies that are spatially distant from engineered *S. epidermidis* are more prone to such toxin outbreaks.

![Fig.34 A Snapshot for Th2 Cytokine in Recurrence 3D Model](media/physicell4/media/image4.png)

Fig.34 A Snapshot for Th2 Cytokine in Recurrence 3D Model

![Fig.35 A Snapshot for Toxin in Recurrence 3D Model](media/physicell4/media/image5.png)

Fig.35 A Snapshot for Toxin in Recurrence 3D Model

Interestingly, **VHH** secretion also manifests as **localized bursts** from *S. epidermidis*, temporally and spatially **responding to nearby *S. aureus* toxin eruptions** — vividly demonstrating the switch’s responsive behavior.

![Fig.36 VHH in Recurrence 3D. Draft still image; published simulation animation(https://static.igem.wiki/teams/5615/wiki/model/3d-lesion-vhh.webp).](media/physicell4/media/image6.png)

Fig.36 VHH in Recurrence 3D. Draft still image; [published simulation animation](https://static.igem.wiki/teams/5615/wiki/model/3d-lesion-vhh.webp).

[Fig.37: VHH snapshot in the Recurrence 3D model](https://static.igem.wiki/teams/5615/wiki/model/snapshot3.webp).

Note: The S.aureus colonies at the two right-side outbreak sites are mixed with S.epidermidis

### Results- Optimal Regions for the Functional Contributions 

When we adjusted the parameters of the *Recurrence* model, we found that **changes in the population sizes** of *S.aureus* and the engineered *S.epidermis* were **highly sensitive to** the **proportion of functional modules** in the **engineered strain**.

Reviewing the recurrence rule set shows that ***S. aureus*** secretes a **toxin** that **stimulates cytokines**, thereby **reducing AMP** production locally and creating **a more favorable niche**. Conversely, the **engineered strain** secretes **VHH** to **reduce** inflammatory **cytokines**, which indirectly **promotes** local **AMP** production that is **detrimental to the engineered strain itself**.

Therefore, if **VHH** secretion or action is **too strong**, the engineered S. epidermidis strain can be rapidly **killed by the resulting AMP surge**. Conversely, if VHH secretion or action is **too weak**, the strain **cannot achieve a therapeutic effect**. Similarly, **AIP-C secretion requires precise adjustment**: increasing AIP-C reduces S. aureus invasion and toxin secretion (which is beneficial for the survival of the engineered strain), but excessive AIP-C can reduce the effectiveness of AIP_sa and limit its ability to promote AMP production.

![Fig.38 Pathways for S.aureus to Occupy Advantageous Ecological Niches](media/physicell4/media/image7.png)

Fig.38 Pathways for *S.aureus* to Occupy Advantageous Ecological Niches

![Fig.39 Pathways for engineered S.epidermis to be introduced to unfavorable ecological niches](media/physicell4/media/image8.png)

Fig.39 Pathways for engineered *S.epidermis* to be introduced to unfavorable ecological niches

Therefore, **both VHH and AIP-C require fine-tuning to find an optimal parameter regime.**

After preliminary calibration, we performed a uniform scan over **three key parameters**:

| <span class="mark">Parameter</span> | <span class="mark">Range</span> | <span class="mark">Number of samples</span> |
|:--:|:--:|:--:|
| <span class="mark">VHH_to_cytokine_degradation_rate</span> | <span class="mark">\[0.1,1\]</span> | <span class="mark">10</span> |
| <span class="mark">Saturation value for the increase of AIP_sa on AIP-C secretion for engr\_*S.epidermis*</span> | <span class="mark">\[0.5,3\]</span> | <span class="mark">10</span> |
| <span class="mark">Saturation value for the increase of AIP_sa on VHH secretion for engr\_*S.epidermis*</span> | <span class="mark">\[0.01,0.5\]</span> | <span class="mark">10</span> |

Next, we built a 10×10×10 parameter space. To map the optimal surface, we defined a **scoring metric** for each simulation. The ideal target state is as follows: **S. aureus** **(A)** s**hould remain constant** (Ai ≈ Ao) with only minor fluctuations; sustained increases or large spikes are undesirable. **Engineered S. epidermidis** **(E) should not decrease** — it should preferably remain constant or increase slightly (Ei ≥ Eo); sustained decline or large drops are undesirable.

Let the time series have $N$ points, $t_0 < t_1 < \cdots < t_{N-1}$. Define $A_i$ as the *S. aureus* count and $E_i$ as the engineered *S. epidermidis* count at $t_i$; the total span is $T = t_{N-1} - t_0$. To avoid division by zero, set $\varepsilon = 10^{-9}$.

We combined three components into a composite penalty score to measure deviation from the ideal:

**1. Average Baseline Deviation** (only penalize undesired direction)

For A (only upward deviations):

```math
D_{A} = \frac{1}{N} \sum_{i=0}^{N-1} \max \left( 0, \frac{A_{i} - A_{0}}{A_{0} + \varepsilon} \right)
```

For E (only downward deviations):

```math
D_{E} = \frac{1}{N} \sum_{i=0}^{N-1} \max \left(0, \frac{E_0 - E_i}{E_0 + \varepsilon} \right)
```

**2. Net trend** (endpoint vs. start; captures overall direction)

For A (penalize final \> initial):

```math
T_{A} = \max \left( 0, \frac{A_{N-1} - A_{0}}{A_{0} + \varepsilon} \right)
```

For E (penalize final \< initial):

```math
T_{E} = \max \left(0, \frac{E_0 - E_{N-1}}{E_0 + \varepsilon} \right)
```

**3. Cumulative single-step direction violations** (local directionality)

For A (sum of all positive single-step increases):

```math
V_{A} = \frac{\sum_{i=1}^{N-1} \max(0, A_i - A_{i-1})}{A_0(N - 1) + \varepsilon}
```

For E (sum of all negative single-step decreases):

```math
V_{E} = \frac{\sum_{i=1}^{N-1} \max(0, E_{i-1} - E_i)}{E_0(N - 1) + \varepsilon}
```

**Overall penalty:**

```math
P = w_{A}^{D} D_{A} + w_{A}^{T} T_{A} + w_{A}^{V} V_{A} + w_{E}^{D} D_{E} + w_{E}^{T} T_{E} + w_{E}^{V} V_{E}
```

We set trend (T) and average deviation (D) with similar importance and direction violations (V) slightly lower:

```math
w_A^D = 1, w_A^T = 1, w_A^V = 0.5
```

```math
w_E^D = 1, w_E^T = 1, w_E^V = 0.5
```

Scores were mapped to \[0,1\] by an inverse transform:

```math
\mathrm{Score} = \frac{1}{1 + P'}
```

![Fig.40 Parameter Scan Score Space](media/physicell4/media/image9.png)

Fig.40 Parameter Scan Score Space

As shown, we obtained a scored parameter space and highlighted the top 5% scoring points.

![Fig.41 Top View of the Parameter Scan Score Space](media/physicell4/media/image10.png)

Fig.41 Top View of the Parameter Scan Score Space

From the top view of the parameter space, **the optimal points are concentrated where the VHH-to-cytokine-degradation-rate and VHH-secretion-saturation are balanced**, while the effect of AIP-C-secretion-saturation is relatively smaller. A much larger score drop is caused by excessive VHH than by VHH insufficiency, indicating that **overly strong antibody activity perturbs the skin microbiome more severely and may induce immune drift**. (de Wijs et al., 2020).

Consequently, these results guide **downstream strain optimization**: it is critical to **measure VHH_to_cytokine_degradation_rate accurately** and to **directionally tune VHH_secretion_saturation** (for example via P3 promoter engineering or VHH amino-acid sequence modification) so the engineered strain can achieve a suitable compromise between colonization and suppression functions.

de Wijs, L. E. M., Nguyen, N. T., Kunkeler, A. C. M., Nijsten, T., Damman, J., & Hijnen, D. J. (2020). Clinical and histopathological characterization of paradoxical head and neck erythema in patients with atopic dermatitis treated with dupilumab: a case series. *Br J Dermatol*, *183*(4), 745-749. <https://doi.org/10.1111/bjd.18730>
