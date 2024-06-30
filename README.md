###### Paper published on arXiv
### Formalizing the causal interpretation in accelerated failure time models with unmeasured heterogeneity
Authors: **M. Brathovde**, H. Putter, M. Valberg and R.A.J. Post (mari.brathovde (a) medisin.uio.no)

&nbsp;[![License](images/badge-mit.svg)](https://opensource.org/license/mit) [![DOI](images/badge-doi.svg)](https://doi.org/10.48550/arXiv.2409.01983) [![arXiv](images/badge-arxiv.svg)](https://arxiv.org/abs/2409.01983)

In the presence of unmeasured heterogeneity, the hazard ratio for exposure has a complex causal interpretation. To address this, accelerated failure time (AFT) models, which assess the effect on the survival time ratio scale, are often suggested as a better alternative. AFT models also allow for straightforward confounder adjustment. In this work, we formalize the causal interpretation of the acceleration factor in AFT models using structural causal models and data under independent censoring. We prove that the acceleration factor is a valid causal effect measure, even in the presence of frailty and treatment effect heterogeneity. Through simulations, we show that the acceleration factor better captures the causal effect than the hazard ratio when both AFT and proportional hazards models apply. Additionally, we extend the interpretation to systems with time-dependent acceleration factors, revealing the challenge of distinguishing between a time-varying homogeneous effect and unmeasured heterogeneity. While the causal interpretation of acceleration factors is promising, we caution practitioners about potential challenges in estimating these factors in the presence of effect heterogeneity.

This repository reproduces the results presented in the paper.

## Citation
If you find the paper or code useful in your work, please cite as BibTeX
```bibtex
@misc{brathovde2024formalizingcausalinterpretationaccelerated,
      title={Formalizing the causal interpretation in accelerated failure time models with unmeasured heterogeneity}, 
      author={Mari Brathovde and Hein Putter and Morten Valberg and Richard A. J. Post},
      year={2024},
      eprint={2409.01983},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2409.01983}, 
}
```
