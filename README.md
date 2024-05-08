# Groundhog Day project

Project Universally Unique Identifier: 9811beb3-2dc5-481f-8b0f-2e8c12936507

**Repository:** `groundhog_day_targets`

In this project, the R scripts are used within the `targets` workflow.

The purpose of the project is to examine the mood drift described by Jangraw et al. (2023) in the context of the Autobiographical Probabilistic Reversal Learning task. The hypothesis is that, when participants are more ingaged with the task, as it happens in the APRL task, then the mood drit disappears (or reverses).

**TODO**

1. Evalutate whether there is a increase of mood during a session.
2. Test whether the APRL impact the post-mood level.
3. Formulate and test an instant mood model, which include an RPE component. MLE version based on a single session.
4. Verify whether is possible formulate a Bayesian hierarchical model of momentary mood, which takes into consideration all subjects and sessions.

Wed Nov  1 08:03:30 CET 202

The current targets workflow generates the data frame `params_happiness_df`, which includes the estimates of the parameters of the model for each subject and session. I need to remove the outliers from this data frame.
