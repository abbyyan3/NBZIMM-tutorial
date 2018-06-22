---
title: home
---

# NBZIMM: Negative Binomial and Zero-Inflated Mixed Models

This R package provides functions for setting up and fitting negative binomial mixed models and zero-inflated negative binomial, Gaussian and Poisson models. These functions allow for mutiple and correlated group-specific (random) effects and various types of within-group correlation structures as described in the core package nlme, and return objects that can be summarized by functions in nlme. The methods can be used to analyze overdispersed and zero-inflated count or continuous responses with multilevel data structures (for example, clustered and longitudinal studies).  


Author: Nengjun Yi nyi@uab.edu

Maintainer: Nengjun Yi nyi@uab.edu

<div class="toc" markdown="1">
## The following links are some example codes used for manuscripts with NBZIMMs:

{% for lesson in site.pages %}
{% if lesson.nav == true %}- [{{ lesson.title }}]({{ lesson.url | absolute_url }}){% endif %}
{% endfor %}
</div>



**License**: GPL

**Author**: Nengjun Yi <nyi@uab.edu>

**Maintainer**: Nengjun Yi <nyi@uab.edu>; Xinyan Zhang <xzhang@georgiasouthern.edu>


