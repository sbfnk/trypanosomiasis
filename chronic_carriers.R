library('data.table')
library('XLConnect')

wb <- loadWorkbook("~/Research/Analysis/Sleeping_Sickness/Chronic/Village input data.xlsx")
village_screening <- data.table(readWorksheet(wb, "Sheet1"))

village_dirs <- list.files(path = path.expand("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/"))
village_ids <- as.integer(gsub("^village ", "", village_dirs))
village_ids <- village_ids[order(village_ids)]

village_cases <- data.table(village.number = integer(0), month = integer(0), cases = integer(0), stage = integer(0))
for (id in village_ids)
{
    for (stage in c(1, 2))
    {
        temp <- data.table(read.csv(paste("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/village ",
                                          id, "/ps", stage, "_", id, ".csv", sep = ""), header = F))
        setnames(temp, 1:2, c("month", "cases"))
        village_cases <- rbind(village_cases, data.table(village.number = rep(id, nrow(temp)),
                                                         month = temp[, month],
                                                         cases = temp[, cases],
                                                         stage = rep(stage, nrow(temp))))
    }
}



sim_chronic <- function(tmax,  params)
{
    S <- 
}

# Differential Equations
next susceptible[1..b] =susceptible[i] - S_I1[i] - S_I1c[i] +  I1c_S[i] + ps1_removal + ps2_removal
next stage1[1..b] = stage1[i] + S_I1[i] - I1_I2[i] - ps1_removal
next stage1c[1..b] = stage1c[i] + S_I1c[i] - I1c_S[i]
next stage2[1..b] = stage2[i] + I1_I2[i] - I2_D[i] - ps2_removal
next deaths[1..b] = deaths[i] + I2_D[i]

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Initial conditions}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
INIT susceptible[1..b] = S_0
INIT stage1[1..b] = i1_0[i]
INIT stage1c[1..b] = i1c_0[i]
INIT stage2[1..b] = i2_0[i]
INIT deaths[1..b] = d_0

{Dimensions of matrix}
b = 2970						; rows = batchruns with different parameter values

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Model starting conditions}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
S_0 = N
; initial susceptibles, assuming 

i1_0[1..b] = actual1[i] - observed1[i]
; initial stage1 prevalence in model after first as

i1c_0[1..b] = actual1c[i] - observed1c[i]
; initial stage1c prevalence in model after first as

i2_0[1..b] = actual2[i] - detected2_1
; initial stage2 prevalence in model after first as

d_0 = 0
; initial deaths

{actual prevalent cases}
{--------------------------------------------------------------------------------------------------------------------------------}
actual1[1..b] = POISSON(observed1[i] / sigma_start[i])			; estimated actual number of stage1 cases
actual1c[1..b] = POISSON(observed1c[i] / (sigma_start[i]*alpha[i]))	; estimated actual number of stage1c cases
actual2[1..b] = POISSON(detected2_1 / sigma_start[i])			; estimated actual number of stage2 cases

sigma_start[1..b] = RANDOM(min(1.05,1),1)					; probability of observing cases at start
sigma_end[1..b] = RANDOM(min(0.72,1),1)					; probability of observing cases at end

{initial observed cases}
{--------------------------------------------------------------------------------------------------------------------------------}
detected1_1 = 1													; number of stage 1 (chronic and pathogenic) detected during as
observed1[1..b] = if STOCH=0 then (detected1_1*proportion_1[i]) else binomial(proportion_1[i], detected1_1)		; number of stage 1 pathogenic detected during as
observed1c[1..b] = if STOCH=0 then (detected1_1*proportion_1c[i]) else binomial(proportion_1c[i], detected1_1)	; number of stage 1 chronic detected during as
detected2_1 = 0													; number of stage 2 detected during as

proportion_1c[1..b] = I1c_e[i]/(I1_e[i] + I1c_e[i])
proportion_1[1..b] = I1_e[i]/(I1_e[i] + I1c_e[i])

{equilibrium prevalence}
{--------------------------------------------------------------------------------------------------------------------------------}
I1_e[1..b] = (lambda[i]*(1-pc[i])*S_0) / (r1)			; equilibrium stage 1
I1c_e[1..b] = (lambda[i]*pc[i]*S_0) / (gamma)		; equilibrium stage 1c
I2_e[1..b] = (r1*I1_e[i]) / (r2)				; equilibrium stage 2

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Infection related parameters - in monthly units}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{rates of disease progression}
r1 = 1 - exp(-0.0019*30.42)				; flow from stage 1 - stage 2
r2 = 1 - exp(-0.0040*30.42)				; flow from stage 2 - deaths (due to HAT)

alpha[1..b] = #alpha(i)					; relative detectability of chronic cf stage1
lambda[1..b] = #lambda(i)				; incident rate of infection
pc[1..b] = #pc(i)						; proportion of infections that are chronic

foi_1c[1..b] = lambda[i]*pc[i]				; force of infection to chronic
foi_1[1..b] = lambda[i]*(1 - pc[i])				; force of infection to stage1

gamma = 1 - exp(-1/120)					; rate at which chronic infection resolves

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Demographic parameters}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
N = 1894						; village specific estimated population at first active screening

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Transitions}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
STOCH=1	; STOCH=0 means deterministic, STOCH=1 is stochastic

S_I1[1..b] = if STOCH=0 then (foi_1[i]*DT*susceptible[i]) else binomial(foi_1[i]*DT,susceptible[i])
;transition from susceptible to pathogenic stage 1 infection

S_I1c[1..b] = if STOCH=0 then (foi_1c[i]*DT*susceptible[i]) else binomial(foi_1c[i]*DT,susceptible[i])
; transition from susceptible to chronic stage 1 infection

I1c_S[1..b] = if STOCH=0 then (gamma*DT*stage1c[i]) else binomial(gamma*DT, stage1c[i])
; transition from chronic stage 1 to susceptible infection

I1_I2[1..b] = if STOCH=0 then (r1*DT*stage1[i]) else (if stage1[i]>0 then binomial(r1*DT,stage1[i]) else 0)
; transition from pathogenic stage 1 to stage 2

I2_D[1..b] = if STOCH=0 then (r2*DT*stage2[i]) else (if stage2[i]>0 then binomial(r2*DT,stage2[i]) else 0)
; transition from stage 2 to HAT specific death

{removal of pathogenic stage 1 and stage 2 due to passive screening}
{--------------------------------------------------------------------------------------------------------------------------------}
ps1 = #ps1_1(time)					; read in ps1 dataset - passive detection of stage 1 cases
ps2 = #ps2_1(time)					; read in ps2 dataset - passive detection of stage 2 cases

ps1_removal = ps1					; ps1 removal of cases passively detected
ps2_removal = ps2					; ps2 removal of cases passively detected

{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{Analysis}
{----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------}
{predicted end prevalence}
{--------------------------------------------------------------------------------------------------------------------------------}
detected1_2 = 0						; observed end prevalence of stage 1 (chronic plus pathogenic) cases

predicted_end1[1..b] = if STOCH=0 then (if time=STOPTIME then ((stage1[i]*sigma_end[i]) + (stage1c[i]*sigma_end[i]*alpha[i])) else 0) else (if time=STOPTIME then (binomial(sigma_end[i]*DT,stage1[i]) + binomial(sigma_end[i]*alpha[i]*DT,stage1c[i])) else 0)
							; predicted end prevalence of stage1 (chronic plus pathogenic) cases

{sum of squares}
{--------------------------------------------------------------------------------------------------------------------------------}
sum_squares_1[1..b] = if neg_I1[i]=1000000 then (if time=STOPTIME then 1000000 else 0) else (if time=STOPTIME then ((predicted_end1[i] - detected1_2 )^2) else 0)  
{ FC PhD OPTION: sum_squares_1[1..b] = if neg_I1[i]=1000000 then (if time=STOPTIME then 0 else 0) else (if time=STOPTIME then (if predicted_end1[i] = detected1_2 then 1 else 0) else 0)  }

{parameter plot dummy variable}
{--------------------------------------------------------------------------------------------------------------------------------}
run = 1

{penalise negative cases}
{--------------------------------------------------------------------------------------------------------------------------------}
neg_I1[1..b] = if neg_I1[i]=0 and stage1[i]<0 then 1000000 else neg_I1[i]

