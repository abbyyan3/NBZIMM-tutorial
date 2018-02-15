---
title: home
---

# BhGLM: Bayesian hierarchical GLMs and survival models, with application to Genetics and Epidemiology 

This package provides functions for setting up and fitting various Bayesian hierarchical models (generalized linear models (GLMs), Cox survival models, negative binomial models, zero-inflated negative binomial models, and ordered logistic or probit regressions), for numerically and graphically summarizing the fitted models, and for evaluating the predictive performance. Several types of priors on the coefficients can be used: double-exponential, Student-t, mixture double-exponential, and mixture t. The methods can be used to analyze not only general data but also large-scale molecular data (i.e., detecting disease-associated genes or variants, predictive and prognostic modeling of diseases and traits, and microbiome data analysis).
       
The BhGLM package is actively and openly developed on GitHub: 

(https://github.com/abbyyan3/BhGLM/)

<div class="toc" markdown="1">
## The following links are some example codes used for manuscripts with BhGLM:

{% for lesson in site.pages %}
{% if lesson.nav == true %}- [{{ lesson.title }}]({{ lesson.url | absolute_url }}){% endif %}
{% endfor %}
</div>

Author: Nengjun Yi <nyi@uab.edu>
Maintainer: Nengjun Yi <nyi@uab.edu>; Xinyan Zhang <xzhang@georgiasouthern.edu>
License: GPL

