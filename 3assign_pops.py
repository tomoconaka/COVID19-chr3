import sys
import subprocess
import pkg_resources
import argparse

required = {'sklearn', 'gcsfs>=0.2.2'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing])

import pandas as pd
import numpy as np
from collections import defaultdict, namedtuple, OrderedDict
from typing import *
import random
from sklearn.ensemble import RandomForestClassifier

def load_merge_data(refscores, ref_info, datascores, outdir):
    ref = pd.read_table(refscores, header=0, sep='\t', compression='gzip')
    data = pd.read_table(datascores, header=0, sep='\t')
    ref_info = pd.read_table(ref_info, header=0, sep=';')
    ref_info = ref_info[['Sample', 'Population']]

    d = {'CHB': 'EAS',
         'JPT': 'EAS',
         'CHS': 'EAS',
         'CDX': 'EAS',
         'KHV': 'EAS',
         'CEU': 'EUR',
         'TSI': 'EUR',
         'FIN': 'EUR',
         'GBR': 'EUR',
         'IBS': 'EUR',
         'YRI': 'AFR',
         'LWK': 'AFR',
         'GWD': 'AFR',
         'MSL': 'AFR',
         'ESN': 'AFR',
         'ASW': 'AFR',
         'ACB': 'AFR',
         'MXL': 'AMR',
         'PUR': 'AMR',
         'CLM': 'AMR',
         'PEL': 'AMR',
         'GIH': 'SAS',
         'PJL': 'SAS',
         'BEB': 'SAS',
         'STU': 'SAS',
         'ITU': 'SAS'}

    ref_info['SuperPop'] = ref_info['Population'].map(d)

    ref_info.to_csv(outdir + 'G1000_pops.txt', sep='\t', index=False)

    print("ref columns: " + str(ref.columns.values))
    print("data columns: " + str(data.columns.values))
    print("ref_info columns: " + str(ref_info.columns.values))

    ref_merge = pd.merge(left=ref, right=ref_info, left_on='s', right_on='Sample', how='inner')
    print("ref_merge columns: " + str(ref_merge.columns.values))

    data_ref = pd.concat([ref_merge, data])

    return data_ref


def assign_population_pcs(
        pop_pc_pd: pd.DataFrame,
        num_pcs: int,
        known_col: str = 'SuperPop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth'
) -> Tuple[pd.DataFrame, RandomForestClassifier]:
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.
    :param Table pop_pc_pd: Pandas dataframe containing population PCs as well as a column with population labels
    :param str known_col: Column storing the known population labels
    :param str pcs_col: Columns storing the PCs
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int num_pcs: number of population PCs on which to train the model
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Dataframe containing sample IDs and imputed population labels, trained random forest model
    :rtype: DataFrame, RandomForestClassifier
    """

    # Expand PC column
    # pop_pc_pd = expand_pd_array_col(pop_pc_pd, pcs_col, num_pcs, 'PC')
    pc_cols = ['PC{}'.format(i + 1) for i in range(num_pcs)]
    print(pc_cols)
    # pop_pc_pd[pc_cols] = pd.DataFrame(pop_pc_pd[pcs_col].values.tolist())[list(range(num_pcs))]
    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    N = len(train_data)
    print(train_data.shape)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].as_matrix()
        training_set_pcs = train_fit[pc_cols].as_matrix()
        evaluation_set_pcs = evaluate_fit[pc_cols].as_matrix()

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    print('Classifying data')
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].as_matrix())
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].as_matrix())
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    print('probs shape ' + str(probs.shape))
    print('pop_pc_pd shape ' + str(pop_pc_pd.shape))
    print(probs.iloc[:3,])
    print(pop_pc_pd.iloc[:3,])
    pop_pc_pd = pd.concat([pop_pc_pd.reset_index(drop=True), probs.reset_index(drop=True)], axis=1)
    print(pop_pc_pd.shape)
    print(pop_pc_pd.iloc[:3, ])
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    #pop_pc_pd = pop_pc_pd.drop(pc_cols, axis='columns')
    print(pop_pc_pd.shape)

    return pop_pc_pd, pop_clf


def main(args):
    # error handling -----------------------------------------------------------
    if not args.ref_scores:
        raise Exception("Reference scores file not selected")

    if not args.data_scores:
        raise Exception("Data scores file not selected")

    if not args.out_dir:
        raise Exception("Output directory where files will be saved is not specified")

    data_ref = load_merge_data(args.ref_scores, args.ref_info, args.data_scores, args.out_dir)

    rf_thresh = [0.5, 0.8]
    for threshold in rf_thresh:
        pcs_df, clf = assign_population_pcs(pop_pc_pd=data_ref, num_pcs=20, min_prob=threshold)

        data_pops = pcs_df.loc[pcs_df['SuperPop'].isnull()]
        data_pops['pop'].value_counts()
        cols = ['s', 'pop'] + [f'prob_{i}' for i in ["EUR", "EAS", "AMR", "AFR", "SAS"]] + [f'PC{i}' for i in range(1, 21)]
        data_pops_df = data_pops[cols]
        print(data_pops_df)

        data_pops_df.to_csv('{}cases_controls_pca_sup_pops_{}_probs.txt'.format(args.out_dir, threshold),
                            sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_scores', type=str, required=True)
    parser.add_argument('--ref_info', default='gs://dsge-covid19-data/G1000/20130606_sample_info.csv')
    parser.add_argument('--data_scores', type=str, required=True)
    parser.add_argument('--out_dir', type=str, required=True)

    args = parser.parse_args()
    main(args)


# EXAMPLE
# hailctl dataproc submit hail assign_pops.py
# --ref_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/1000G_scores.txt.bgz
# --data_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/COV_ILLUMINA_15102020.chr0.pos0.removed.qc_scores.tsv
# --out_dir gs://dsge-covid19-data/COV_ILLUMINA_15102020/pca/