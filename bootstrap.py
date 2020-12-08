import numpy as np
import pandas as pd


D_samples = pd.read_csv("data/samples_D.tsv.gz", sep="\t", index_col=[0])
R_samples = pd.read_csv("data/samples_R.tsv.gz", sep="\t", index_col=[0])


def bootstrap_per_group_count(row, cohort="D"):
    tweet_count = pd.read_csv("results/{}_number_of_tweets.tsv".format(cohort), sep="\t", index_col="user_id")
    return row.apply(lambda x: tweet_count.loc[x, "count"]).sum()


bootstrap = pd.DataFrame(index=pd.Int64Index(np.arange(1, 10001), name="num_runs"),
                         columns=["$D$", "$R$"])
bootstrap.loc[:, "$D$"] = D_samples.apply(bootstrap_per_group_count, cohort="D", axis=1)
bootstrap.loc[:, "$R$"] = R_samples.apply(bootstrap_per_group_count, cohort="R", axis=1)
bootstrap.to_csv("bootstrap/tweets_per_run.tsv", sep="\t")


def bootstrap_per_group(row, cohort="D", without_FPP=False):
    fpp = "_without_FPP" if without_FPP else ""
    fn = "results/{}_number_of_tweets_with_CDS{}.tsv".format(cohort, fpp)
    CDS_count = pd.read_csv(fn, sep="\t", index_col="user_id")
    tweet_count = pd.read_csv("results/{}_number_of_tweets.tsv".format(cohort), sep="\t", index_col="user_id")
    return row.apply(lambda x: CDS_count.loc[x, "count"]).sum() / row.apply(lambda x: tweet_count.loc[x, "count"]).sum()


D_P = D_samples.apply(bootstrap_per_group, cohort="D", axis=1)
R_P = R_samples.apply(bootstrap_per_group, cohort="R", axis=1)
RP = D_P / R_P
RP.to_csv("bootstrap/relative_prevalence.tsv", sep="\t")
del D_P, R_P, RP

D_P_nonFPP = D_samples.apply(bootstrap_per_group, cohort="D", without_FPP=True, axis=1)
R_P_nonFPP = R_samples.apply(bootstrap_per_group, cohort="R", without_FPP=True, axis=1)
RP_nonFPP = D_P_nonFPP / R_P_nonFPP
RP_nonFPP.to_csv("bootstrap/relative_prevalence_without_FPP.tsv", sep="\t")
del D_P_nonFPP, R_P_nonFPP, RP_nonFPP

D_samples_male = pd.read_csv("data/samples_D_male.tsv.gz", sep="\t", index_col=[0])
R_samples_male = pd.read_csv("data/samples_R_male.tsv.gz", sep="\t", index_col=[0])
D_P_male = D_samples_male.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_male = R_samples_male.apply(bootstrap_per_group, cohort="R", axis=1)
RP_male = D_P_male / R_P_male
RP_male.to_csv("bootstrap/relative_prevalence_male.tsv", sep="\t")
del D_samples_male, R_samples_male, D_P_male, R_P_male, RP_male

D_samples_female = pd.read_csv("data/samples_D_female.tsv.gz", sep="\t", index_col=[0])
R_samples_female = pd.read_csv("data/samples_R_female.tsv.gz", sep="\t", index_col=[0])
D_P_female = D_samples_female.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_female = R_samples_female.apply(bootstrap_per_group, cohort="R", axis=1)
RP_female = D_P_female / R_P_female
RP_female.to_csv("bootstrap/relative_prevalence_female.tsv", sep="\t")
del D_samples_female, R_samples_female, D_P_female, R_P_female, RP_female

D_samples_teens = pd.read_csv("data/samples_D_teens.tsv.gz", sep="\t", index_col=[0])
R_samples_teens = pd.read_csv("data/samples_R_teens.tsv.gz", sep="\t", index_col=[0])
D_P_teens = D_samples_teens.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_teens = R_samples_teens.apply(bootstrap_per_group, cohort="R", axis=1)
RP_teens = D_P_teens / R_P_teens
RP_teens.to_csv("bootstrap/relative_prevalence_teens.tsv", sep="\t")
del D_samples_teens, R_samples_teens, D_P_teens, R_P_teens, RP_teens

D_samples_twenties = pd.read_csv("data/samples_D_twenties.tsv.gz", sep="\t", index_col=[0])
R_samples_twenties = pd.read_csv("data/samples_R_twenties.tsv.gz", sep="\t", index_col=[0])
D_P_twenties = D_samples_twenties.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_twenties = R_samples_twenties.apply(bootstrap_per_group, cohort="R", axis=1)
RP_twenties = D_P_twenties / R_P_twenties
RP_twenties.to_csv("bootstrap/relative_prevalence_twenties.tsv", sep="\t")
del D_samples_twenties, R_samples_twenties, D_P_twenties, R_P_twenties, RP_twenties

D_samples_thirties = pd.read_csv("data/samples_D_thirties.tsv.gz", sep="\t", index_col=[0])
R_samples_thirties = pd.read_csv("data/samples_R_thirties.tsv.gz", sep="\t", index_col=[0])
D_P_thirties = D_samples_thirties.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_thirties = R_samples_thirties.apply(bootstrap_per_group, cohort="R", axis=1)
RP_thirties = D_P_thirties / R_P_thirties
RP_thirties.to_csv("bootstrap/relative_prevalence_thirties.tsv", sep="\t")
del D_samples_thirties, R_samples_thirties, D_P_thirties, R_P_thirties, RP_thirties

D_samples_adults = pd.read_csv("data/samples_D_adults.tsv.gz", sep="\t", index_col=[0])
R_samples_adults = pd.read_csv("data/samples_R_adults.tsv.gz", sep="\t", index_col=[0])
D_P_adults = D_samples_adults.apply(bootstrap_per_group, cohort="D", axis=1)
R_P_adults = R_samples_adults.apply(bootstrap_per_group, cohort="R", axis=1)
RP_adults = D_P_adults / R_P_adults
RP_adults.to_csv("bootstrap/relative_prevalence_adults.tsv", sep="\t")
del D_samples_adults, R_samples_adults, D_P_adults, R_P_adults, RP_adults


def bootstrap_per_category(row, cohort="D", without_FPP=False):
    fpp = "_without_FPP" if without_FPP else ""
    fn = "results/{}_number_of_tweets_with_CDS{}_per_category.tsv".format(cohort, fpp)
    CDS_count = pd.read_csv(fn, sep="\t", index_col="user_id")
    tweet_count = pd.read_csv("results/{}_number_of_tweets.tsv".format(cohort), sep="\t", index_col="user_id")
    return row.apply(lambda x: CDS_count.loc[x, :]).sum() / row.apply(lambda x: tweet_count.loc[x, "count"]).sum()


D_P_category = D_samples.apply(bootstrap_per_category, cohort="D", axis=1)
R_P_category = R_samples.apply(bootstrap_per_category, cohort="R", axis=1)
RP_category = D_P_category / R_P_category
RP_category.to_csv("bootstrap/relative_prevalence_category.tsv", sep="\t")
del D_P_category, R_P_category, RP_category

D_P_nonFPP_cat = D_samples.apply(bootstrap_per_category, cohort="D", without_FPP=True, axis=1)
R_P_nonFPP_cat = R_samples.apply(bootstrap_per_category, cohort="R", without_FPP=True, axis=1)
RP_nonFPP_cat = D_P_nonFPP_cat / R_P_nonFPP_cat
RP_nonFPP_cat.to_csv("bootstrap/relative_prevalence_category_without_FPP.tsv", sep="\t")
del D_P_nonFPP_cat, R_P_nonFPP_cat, RP_nonFPP_cat


def bootstrap_per_phrase(row, cohort="D"):
    fn = "results/{}_number_of_tweets_with_CDS_per_schema.tsv".format(cohort)
    CDS_count = pd.read_csv(fn, sep="\t", index_col="user_id")
    tweet_count = pd.read_csv("results/{}_number_of_tweets.tsv".format(cohort), sep="\t", index_col="user_id")
    return row.apply(lambda x: CDS_count.loc[x, :]).sum() / row.apply(lambda x: tweet_count.loc[x, "count"]).sum()


D_P_phrase = D_samples.apply(bootstrap_per_phrase, cohort="D", axis=1)
R_P_phrase = R_samples.apply(bootstrap_per_phrase, cohort="R", axis=1)
RP_phrase = D_P_phrase / R_P_phrase
RP_phrase.to_csv("bootstrap/relative_prevalence_phrase.tsv", sep="\t")
del D_P_phrase, R_P_phrase, RP_phrase


def bootstrap_CDS(row, cohort="D"):
    fn = "results/{}_number_of_tweets_with_CDS_per_schema.tsv".format(cohort)
    CDS_count = pd.read_csv(fn, sep="\t", index_col="user_id")
    tweet_count = pd.read_csv("results/{}_number_of_tweets.tsv".format(cohort), sep="\t", index_col="user_id")
    return row.apply(lambda x: CDS_count.loc[:, x]).sum().sum() / tweet_count["count"].sum()


CDS_samples = pd.read_csv("data/samples_CDS.tsv", sep="\t", index_col=[0])
D_P_CDS = CDS_samples.apply(bootstrap_CDS, cohort="D", axis=1)
R_P_CDS = CDS_samples.apply(bootstrap_CDS, cohort="R", axis=1)
RP_CDS = D_P_CDS / R_P_CDS
RP_CDS.to_csv("bootstrap/relative_prevalence_CDS.tsv", sep="\t")
del D_P_CDS, R_P_CDS, RP_CDS

_CDS = pd.read_csv("data/list_of_CDS.tsv", sep="\t", index_col="markers")
for cat in _CDS.categories.unique():
    CDS_samples = pd.read_csv("data/samples_CDS_{}.tsv.gz".format(cat.replace(" ", "-")),
                              sep="\t", index_col="num_runs")
    D_P_CDS = CDS_samples.apply(bootstrap_CDS, cohort="D", axis=1)
    R_P_CDS = CDS_samples.apply(bootstrap_CDS, cohort="R", axis=1)
    RP_CDS = D_P_CDS / R_P_CDS
    RP_CDS.to_csv("bootstrap/relative_prevalence_CDS_{}.tsv".format(cat.replace(" ", "-")),
                  sep="\t")
    del D_P_CDS, R_P_CDS, RP_CDS
