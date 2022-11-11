# HF prediction by Machine Learning

## Preprocessing steps

1. Remove samples:
    i. PrevalentHFAIL == 1; Event == 1; Event_time < 0
    ii. NA's in Event or Event_time
2. Filter taxa:
    i. By counts
    ii. Aglomerating by Species
3. Calculate alpha diversity (7 different scores) from all species
