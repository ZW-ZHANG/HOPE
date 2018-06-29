# HOPE
This is a sample implementation of "[Asymmetric Transitivity Preserving Graph Embedding](http://www.kdd.org/kdd2016/papers/files/rfp0184-ouA.pdf)"(KDD 2016).

### Requirements
```
MATLAB R2014a
``` 

### Usage
run `embed_main.m` with matlab

```
Input:
	A: N*N adjacency matrix (sparse)

	K: dimensionality of embedding space

	beta: decaying constant, default is 0.5 / spectral radius
    
Output:
    U: N*K left embedding matrix

	V: N*K right embedding matrix

	The high-order proximity (katz) matrix is approximated by U * V'

```

### Cite
If you find this code useful, please cite our paper:
```
@inproceedings{ou2016asymmetric,
  title={Asymmetric transitivity preserving graph embedding},
  author={Ou, Mingdong and Cui, Peng and Pei, Jian and Zhang, Ziwei and Zhu, Wenwu},
  booktitle={Proceedings of the 22nd ACM SIGKDD international conference on Knowledge discovery and data mining},
  pages={1105--1114},
  year={2016},
  organization={ACM}
}
```
