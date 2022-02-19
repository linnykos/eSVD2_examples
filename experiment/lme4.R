library(lme4)
# from https://stats.stackexchange.com/questions/122009/extracting-slopes-for-cases-from-a-mixed-effects-model-lme4
data(sleepstudy)
fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
head(sleepstudy)

summary(fm1)

coef(fm1)
coef(fm1)$Subject
coef(summary(fm1))
coef(summary(fm1))[ , "Estimate"]
lme4::ranef(fm1)$Subject
lme4::ranef(fm1)

res <- predict(fm1, sleepstudy)
coef_mat <- coef(fm1)$Subject
res2 <- sapply(1:nrow(sleepstudy), function(i){
  idx <- sleepstudy[i,"Subject"]
  row_idx <- which(rownames(coef_mat) == idx)
  coef_mat[row_idx,1] + coef_mat[row_idx,2]*sleepstudy[i,"Days"]
})
plot(res, res2, asp = T)

###############

fm1 <- lme4::lmer(Reaction ~ Days + (1|Subject), sleepstudy)
coef(fm1)$Subject

anova(fm1)

###################

# from https://stats.oarc.ucla.edu/r/faq/random-coefficient-poisson-models/
dat <- foreign::read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")
dat$cid <- factor(dat$cid)
m1a <- lme4::glmer(awards ~ 1 + (1 | cid),
                   data = dat,
                   family = stats::poisson(link = "log"))
coef(summary(m1a))
m0 <- lme4::glmer(awards ~ -1 + (1 | cid),
                   data = dat,
                   family = stats::poisson(link = "log"))
coef(m0)$cid
mean(coef(m0)$cid[,1])
coef(summary(m0))

m2 <- lme4::glmer(awards ~ 1 + female + (1 | cid),
            data = dat,
            family = stats::poisson(link = "log"))

# see https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html
anova(m2, m1a)
anova(m1a, m2)

lmtest::lrtest(m2, m1a) # but probably not that great...? see https://www.rdocumentation.org/packages/lmtest/versions/0.9-39/topics/lrtest
