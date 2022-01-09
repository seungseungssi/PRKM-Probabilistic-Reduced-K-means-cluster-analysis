이승훈 ( Seunghoon Lee ) , 송주원 ( Juwon Song ). 2021. 확률적 reduced K-means 군집분석. 응용통계연구, 34(6): 905-922


# PRKM(Probabilistic Reduced K means)

Cluster analysis is one of unsupervised learning techniques used for discovering clusters when there is no prior knowledge of group membership. K-means, one of the commonly used cluster analysis techniques, may fail when the number of variables becomes large. In such high-dimensional cases, it is common to perform tandem analysis, K-means cluster analysis after reducing the number of variables using dimension reduction methods. However, there is no guarantee that the reduced dimension reveals the cluster structure properly. Principal component analysis may mask the structure of clusters, especially when there are large variances for variables that are not related to cluster structure. To overcome this, techniques that perform dimension reduction and cluster analysis simultaneously have been suggested. This study proposes probabilistic reduced K-means, the transition of reduced K-means (De Soete and Caroll, 1994) into a probabilistic framework. Simulation shows that the proposed method performs better than tandem clustering or clustering without any dimension reduction. When the number of the variables is larger than the number of samples in each cluster, probabilistic reduced K-means show better formation of clusters than non-probabilistic reduced K-means. In the application to a real data set, it revealed similar or better cluster structure compared to other methods.


## file descriptions

- PRKM.R : Function for performing Probabilistic Reduced K-means (do not work when k=2, q=1)
- PRKM_for_k2.R : Function for proceeding Probabilistic Reduced K-means (also work when k=2, q=1. But because this function is composed of many for statement, which means elementwise computation(not matrix), it may be slower than PRKM.R, although do not have much difference by experience)
- Thesis_simulation_final.zip : R programming for performing simulation. This also contains functions for GMM, but in my thesis, I used MixAll package for performing GMM
- make_simul_data.R : Function for making simulation data(k=4, q=2)
- Summary_N50.html : Visual summary of simulation result when N = 200(50 observation for each four cluster) and nunber of variables = 10, 50 ,100
- Summary_final.html : Visual summary of simulation result in the thesis(number of variable = 10, 50, 100, N = 200, 400)
- 확률적 Reduced K-means 군집분석 : Thesis for my Master's degree
- EM formula notes for PRKM.pdf : Formula derivation for EM algorithm for PRKM
