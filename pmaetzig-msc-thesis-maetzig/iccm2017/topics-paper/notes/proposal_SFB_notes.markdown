# Notes from SFB Proposal (impaired sentence comprehension)

## Relevant Sections:

1) Summary
2) Research rationale
    1) Current state of understanding and preliminary work
        1) Sentence comprehension and computational modelling in aphasia
        2) Assessment and treatment of sentence comprehension deficits in IWAs and healthy controls
    2) ...
3) Project plan
    1) WP1: Developing a detailed computational model of sentence comprehension of IWAs [...]

## Notes on sections

#### 1) Summary

- individual differences in sentence processing
- linguistic (language-dependent) factors: sentence type / construction
- methodological factors: task used to assess performance
- between-subjects factors: age, IWA/unimpaired
- interactions between factors
- systematicity and limits of individual differences (this is an assumption, albeit a sound one)
- predictions of computational model
- causal links between underlying deficit and improvement in specific sentence comprehension task (proposal specific, not for paper)

#### 2) Research rationale / current state of understanding and preliminary work

(only the first subsection is relevant (because the second subsection is references and contributions), so they are collapsed here)

##### 2.1) Sentence comprehension and computational modelling in aphasia

- sentence comprehension in aphasia: performance of IWAs (individuals with aphasia) subject to a high degree of variability
- [dependent on the measurement/assessment] some show chance performance, others either above or below chance
- healthy controls perform variably and error-prone on non-canonical sentences (Ferreira, 2003)
- performance gradually deteriorates with increasing age (Burcherct et al., 2011) 
- IWAs often show deficiencies in sentence-picture-matching (SPM) or object manipulation (OM) tasks
    - the modelling in our paper could be seen as an instanciation of modelling OM tasks ... the decision that the subject
      and the model have to make are the same (recall dependencies)
    - effect more pronounced with non-canonical word order (in our case, OR vs SR) 

**Assumption**: A deficit that affects syntactic processing mechanisms underlies this behavioural, sentence comprehension deficit in IWAs.

Hypotheses on what exactly that underlying deficit could be (see also Caplan et al., 2015, for an overview):

- slowdown in parsing mechanisms (Burkhardt et al., 2003)
- temporal breakdowns, called *intermittent deficiencies* (Caplan et al., 2015)
- lexical integration problems (Choy & Thompson, 2010)
- resource reduction, regarding working memory (...)

**Role of computational modelling**

- the sentence processing model based on ACT-R framework (cf. Anderson et al., 2004) developed by Lewis & Vasishth (2005) 
- ACT-R framework (as the Lewis & Vasishth, 2005 model) carries with it many assumptions; relevant for the case at hand would be working memory (which is commonly hypothesised to have an influence on sentence processing, cf. Just & Carpenter, 1992), ... (**more assumptions here?**)
- the model developed by/in Lewis & Vasishth (2005) accepts the assumptions of the ACT-R framework and allows sentence comprehension to be treated as a more general information processing task
- Patil et al. (2016) extended the Lewis & Vasishth (2005) model to simulate the performance of impaired sentence comprehension of seven IWAs in a SPM task (data from Hanne et al., 2011)
    - this they did via including (**what is this exactly?**) assumptions about processing slowdown and intermittent deficiencies into the model which was originally created to simulate unimpaired processing
    - best match of the model to the averaged IWA data occured when both slowdown and intermittent deficiencies (i.e., the parameters assumed to correspond to these **deficits [not sure on the word here]**) were included at the same time
    - inter-individual variability could be modelled by manipulating the two parameters corresponding to slowdown and intermittent deficiencies along a [severity] dimension **[better: along a quasi-continuous scale]**
    - **result:** IWAs show variable degrees of deficits along two dimensions, whereas controls show a low range of variability and small values (where small is interpreted as 'less of a deficit / less severe')

**Assumption of this work:** It is possible (and interesting, important) to directly investigate impaired processing by damaging a model of unimpaired sentence processing at points that are hypothesised to contribute to the underlying processing deficit.

Recently, Caplan et al. (2015) addressed the sparsity-of-data problem by creating a corpus of sentence processing data of 61 IWAs (and 41 language unimpaired controls), involving 20 syntactic structures of different complexity. 

**Assumption of this work:** Computational modelling provides a theoretical framework that allows us to investigate, in a theoretically motivated manner, the hypotheses as to what the nature of the underlying aphasic processing deficit might be. [The Lewis & Vasishth (2005) model is especially suited to investigate *resource reduction* because of its implementation in the ACT-R framework.]

The purpose of this work is to contribute to research investigating variability in impaired populations, both within and between participants.
