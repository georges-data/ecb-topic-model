## Loading libraries
library(tidyverse)
library(tidytext)
library(quanteda)
library(quanteda.textstats)
library(tm)
library(textmineR)
library(topicmodels)
library(textcat)
library(topicdoc)
library(lubridate)
library(reshape2)

## Loading data
rawdata <- dataset <- read.table(
  file.choose(),
  sep = "|",
  quote = "",
  fill = TRUE,
  header = TRUE,
  encoding = "UTF-8",
  stringsAsFactors = FALSE
)


### Cleaning Data
cleandata <- dataset %>%
  select(-title, -subtitle) %>%
  rowid_to_column("speechid") %>%
  mutate(date = as.Date(date)) %>%
  mutate(contents = str_replace_all(contents, "[[:punct:]]", " ")) %>%
  mutate(contents = str_replace_all(contents, "[[:digit:]]", " ")) %>%
  mutate(contents = str_replace_all(contents, "\\s+", " ")) %>%
  filter(contents != "") %>%
  mutate(language = textcat(contents))

word_count_df <- strsplit(cleandata$contents, split=" ")
cleandata$word_count <- sapply(word_count_df, length)

language_overview <- cleandata %>% 
  count(language)

cleandata_final <- cleandata %>%
  filter(word_count >= 150) %>%
  filter(language == "english")

view(head(cleandata))

dim(cleandata_final)
sum(cleandata_final$word_count)

## Constructing DTM
ecb_corpus <- corpus(cleandata_final,
                     text_field = "contents")

ecb_tokens <- tokens(ecb_corpus,
                     remove_punct = TRUE,
                     remove_symbols = TRUE,
                     remove_number = TRUE,
                     remove_url = TRUE) %>%
  tokens_select(pattern = stop_words$word, selection = "remove") %>%
  tokens_replace(lexicon::hash_lemmas$token, 
                 lexicon::hash_lemmas$lemma, 
                 valuetype = "fixed") %>%
  tokens_tolower() %>%
  tokens_remove("")

ecb_bigrams <- textstat_collocations(ecb_tokens, min_count = 100)
ecb_bigrams <- ecb_bigrams[1:50,]

ecb_tokens_final <- tokens_compound(ecb_tokens, ecb_bigrams)

ecb_dfm <- dfm(ecb_tokens_final, remove_padding = TRUE)

ECB_tf_idf <- TermDocFreq(ecb_dfm)

ECB_tf_idf %>% filter(term_freq <= 50) %>%
  ggplot(aes(x=term_freq))+
  geom_histogram(binwidth = 1)

ECB_tf_idf_lowerbound <- ECB_tf_idf %>% 
  filter(term_freq <= 10)

ecb_dfm_final <- ecb_dfm %>%
  dfm_remove(ECB_tf_idf_lowerbound$term) %>%
  dfm_remove(c("euro", 
               "central_bank", 
               "ecb", 
               "monetary_policy"))

ecb_dfm_final_df <- convert(ecb_dfm_final, to="data.frame") %>%
  rowid_to_column("posterior_id")

## Running LDA
min_topics <- 2
max_topics <- 10

### Topicmodel package

lda_list_1 <- vector(mode = "list", length=max_topics)
coherence_list_1 <- vector(mode = "list", length=max_topics)

for (topic_no in min_topics:max_topics){
  lda_list_1[[topic_no]] <- LDA(ecb_dfm_final, 
                              k=topic_no, 
                              method="Gibbs", 
                              control=list(iter = 2000, 
                                           seed = 1234, 
                                           verbose = 25))
  
  coherence_list_1[[topic_no]] <- mean(topic_coherence(lda_list_1[[topic_no]], 
                                                       ecb_dfm_final))
}

coherence_list_1[[1]] <- 0

coherence_list_1

coherence_df <- data.frame(matrix(unlist(coherence_list_1), 
                                  nrow=length(coherence_list_1), 
                                  byrow=TRUE)) %>%
  rowid_to_column("topics")

colnames(coherence_df) <- c("topics", "mean_coherence")

coherence_plot <- coherence_df %>%
  filter(topics > 1) %>%
  ggplot(aes(x=topics, y=mean_coherence)) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(2, 10, by = 1))

coherence_plot

topic_overview_1 <- tidy(lda_list_1[[4]], matrix= "beta") %>%
       group_by(topic) %>%
       slice_max(beta, n = 10) %>% 
       ungroup() %>%
       arrange(topic, -beta)

topic_overview_plot_1 <- topic_overview_1 %>%
  mutate(term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(beta, term, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  scale_y_reordered() +
  theme_bw() +
  labs(x = "per-term-per-topic probability (beta)")

topic_overview_plot_1

### Topics on Documents
lda_posterior_per_document <- posterior(lda_list_1[[4]], ecb_dfm_final)

lda_documents <- lda_posterior_per_document$topics

lda_documents <- as.data.frame(lda_documents) %>%
  rowid_to_column("posterior_id")

document_vars <- docvars(ecb_dfm_final) %>%
  rowid_to_column("posterior_id")

lda_documents_final <- inner_join(lda_documents, document_vars)

dim(lda_documents)
dim(document_vars)


## Visualisation
lda_overview <- lda_documents_final %>%
  mutate(quarter = quarter(date)) %>%
  mutate(year = year(date)) %>%
  mutate(year_quarter = year+(quarter-1)/4)

colnames(lda_overview)[2:5] <- c("topic_payments", 
                                 "topic_inflation", 
                                 "topic_economy", 
                                 "topic_financialmarkets")

lda_overview_yearly <- lda_overview %>%
  group_by(year_quarter) %>%
  summarize(payments =  weighted.mean(topic_payments, 
                                                       word_count),
            inflation =  weighted.mean(topic_inflation, 
                                                        word_count),
            economic_growth =  weighted.mean(topic_economy, 
                                                      word_count),
            financial_markets =  weighted.mean(topic_financialmarkets, 
                                                               word_count))

lda_overview_yearly_melted <- lda_overview_yearly %>% 
  melt(id.vars="year_quarter")

view(lda_overview_yearly)

lda_overview_plot_line <- lda_overview_yearly_melted %>%
  filter(year_quarter <= 2021.25) %>%
  ggplot(aes(x=year_quarter, y=value, color=variable)) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(1997, 2022, by = 1)) +
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle=90)) +
  ylim(0,0.8) +
  labs(y = "per-document-per-topic probability (gamma)",
       x = "year",
       color = "topic")

lda_overview_plot_line

lda_overview_plot_bar <- lda_overview_yearly_melted %>%
  filter(year_quarter <= 2021.25) %>%
  ggplot(aes(x=year_quarter, y=value, fill=variable)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_x_continuous(breaks = seq(1997, 2022, by = 1)) +
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle=90))+
  labs(y = "per-document-per-topic probability (gamma)",
       x = "year",
       fill = "topic")

lda_overview_plot_bar
  
### TextmineR package
lda_list_textminer <- vector(mode = "list", length=max_topics)
coherence_list_textminer <- vector(mode = "list", length=max_topics)

for (topic_no in min_topics:max_topics){
  lda_list_textminer[[topic_no]] <- FitLdaModel(
    dtm = ecb_dfm_final,
    k = topic_no,
    iterations = 2000,
    optimize_alpha = TRUE,
    CalcLikelihood = TRUE,
    calc_coherence = TRUE,
    cal_r2 = TRUE
  )
}

textminer_dfm <- convert(ecb_dfm_final, to = "topicmodels")
textminer_lda <- FitLdaModel(
  dtm = ecb_dfm_final,
  k = 2,
  iterations = 100,
  optimize_alpha = TRUE,
  CalcLikelihood = TRUE,
  calc_coherence = TRUE,
  cal_r2 = TRUE)


