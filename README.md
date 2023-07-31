# 本稿の目的   
本稿は、**一般化加法モデル**(GAM)の概要を解説し、それをRで実行する方法を学ぶことを目的とする。本稿の内容は[こちら](https://tsubasayamaguchi-jinrui.github.io/GAM_nyuumon/)から読むことができる。  

GAMは一般化線形モデル(GLMに代表される線形なモデルを拡張し、変数間の関係をより柔軟な形で表現できるようにしたものである。そのため、**GLMで仮定されるような単調増加または単調減少の関係だけでなく、非線形な関係を調べることができる**。

霊長類の行動のような複雑なデータでは変数間の関係が非線形になることがしばしばあるため、GAMは多くの研究で用いられている(e.g., Matsumoto 2017; Taniguchi and Matsumoto-Oda 2018; Hongo et al. 2022)。GLMのように線形性を仮定するモデルがデータに当てはまらない場合には、GAMなどの非線形性を許容するモデルを使用する必要性が生じてくるだろう。

本稿は、Alain Zuurが執筆した”A beginner’s guide to generalized additive models with R”(Zuur 2012)の内容を基に執筆している。本書はなるべく数学的な説明を省きつつ、実際の生態学のデータを用いてGAMについてわかりやすく解説したもので、GAMの入門として非常によい書籍である。より詳細な情報を知りたい場合は原著にアクセスしていただきたい。   

その他に参考にしたのは以下の本である。

- Zuur (2009) Mixed effects models and extensions in ecology with R.  
- James et al. (2013) An Introduction to Statistical Learning with Applications in R.  
- 竹澤 (2009) Rによるノンパラメトリック回帰の入門講義  


**References**    
Hongo S, Nakashima Y, Akomo-Okoue EF, Mindonga-Nguelet FL (2022) Seasonality in daily movement patterns of mandrills revealed by combining direct tracking and camera traps. J Mammal 103:159–168  
James G, D W, Hastie T, Tibshirani R (2013) An Introduction to Statistical Learning with Applications in R. Springer  
Matsumoto T (2017) Developmental changes in feeding behaviors of infant chimpanzees at mahale, tanzania: Implications for nutritional independence long before cessation of nipple contact. Am J Phys Anthropol 163:356–366  
Taniguchi H, Matsumoto-Oda A (2018) Wound healing in wild male baboons: Estimating healing time from wound size. PLoS One 13:  
Zuur AF (2012) A beginner’s guide to generalized additive models with R. Highland Statistics, Newburgh, Scotland  
Zuur AF (2009) Mixed effects models and extensions in ecology with R. Springer, New York, NY  
竹澤邦夫 (2009) Rによるノンパラメトリック回帰の入門講義. メタ・ブレーン  