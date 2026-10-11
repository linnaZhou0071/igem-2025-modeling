# Overview — original wiki draft

> Converted from [Overview.docx](../originals/Overview.docx). Original draft wording, figures and equations are preserved; see [editorial notes](../EDITORIAL_NOTES.md) before interpreting its claims.

Highlights

ODE:

- *In vitro* validation of the system’s stability and therapeutic effect

- Estimation of the final concentration of scFv is 42 mM

- Suggesting replacing scFv with VHH molecules to enhance stability

2D & 3D PhysiCell

- Agent-based cell grammar for microbial interactions and biofilm formation

- Validation of occupation of the epidermal niche and long-term colonization

- Identification of the optimal functional zone for engineered strain

- Overview

Atopic Dermatitis (AD) is a chronic inflammatory skin disorder that can cause severe physical and mental health problems. Our project aims to **develop an engineered strain of *S.epidermis* that can colonize patients’ skin to assist in both treatment and prevention of the recurrence of AD**. A **quorum-sensing responsive switch** enables our engineered strain to interfere with the quorum-sensing system of *S. aureus* to reduce its pathogenicity, as well as secrete IL-4Rα, the Variable domain of Heavy chain of Heavy-chain-only antibody (VHH) to reduce inflammation at the same time.

**Our modeling framework follows a gradual progression toward biological realism.** At the most fundamental level, we conducted **gene circuit stimulation** to illustrate the molecular interactions between engineered *S. epidermidis* and *S. aureus*, focusing on the **theoretical feasibility** of activating the quorum-sensing responsive switch, VHH secreting and disrupting the Agr system through quorum sensing.

Next, we developed a **2D spatial model** to **simulate microbial competition and colonization under different severities of AD**. Since safety concerns prevent us from performing experiments on real animal skin surfaces, this model compensates for the experimental limitation. In addition, we performed a **parameter sweep** to **identify the optimal regions of functional distribution for the engineered strain**. This guided the practical design and further improvement.

Finally, we constructed a **3D model** that goes a step further by **incorporating host skin epithelial cells and T cell immune responses** alongside microbial dynamics. This level most closely resembles the real skin microenvironment, enabling us to assess the long-term control of inflammation and relapse prevention.

As the main highlight, both our 2D and 3D simulations present **a novel approach to modeling microbial communities and their interactions with host cells.** Traditional approaches often rely on population-level equations or bulk dynamics, which lack a clear representation of individual species and their spatial interactions. In contrast, we innovatively adapted a **cell grammar framework** to **define distinct rules for each type of microorganism** and to simulate their competition and coexistence in 2D and 3D spatial environments, while integrating interactions with skin epithelial cells and T cell immune responses. This cross-domain strategy underscores the originality of our modeling framework and provides a more comprehensive understanding of the dynamics underlying our therapeutic design.
