################################################################################
# Title: "data_preparation_descr_analysis.R"                                   #
# Description: Data preparation for the descriptive analysis of data           #
#              collected in the MCID_3 study.                                  #
################################################################################

## Dates of saliva sample collection #------------------------------------------
test.dates <- unique(sal_res$date)[order(as.Date(unique(sal_res$date)))]
time_id <- 1:18
names(time_id) <- test.dates

## Saliva sample test results #-------------------------------------------------
viral_loads <- data.frame(sal_res %>%
                            # Transformation of variables
                            mutate(student_id = as.factor(student_id),
                                   date = as.Date(date),
                                   class = as.factor(class),
                                   weekday = as.factor(weekday),
                                   pathogen = factor(pathogen, levels = pathgs),
                                   panel = as.factor(panel),
                                   ct = as.numeric(ct),
                                   time_id = time_id[match(date, names(time_id))],
                                   test_id = paste0(time_id,"-",substr(student_id,1,1),
                                                    "-", substr(student_id, 2, 5))) %>%
                            # Inclusion of social survey data assessed at baseline
                            left_join(soc_res %>%
                                        filter(time == "Baseline") %>%
                                        dplyr::select(student_id, sex),
                                      join_by(student_id)) %>%
                            rename("gender" = sex) %>%
                            mutate(gender = as.factor(gender)) %>%
                            # If respiratory virus was detected on multiple panels,
                            # consider lowest circular treshold (Ct) value
                            group_by(student_id, date, pathogen) %>%
                            filter(ct == min(ct)) %>%
                            filter(row_number() == 1) %>%
                            ungroup())

## School absenteeism data #----------------------------------------------------
# Absences due to sickness with respiratory symptoms
absence_viral <- abs_sus %>%
  filter(reason == "sickness", symp_respir == 1) %>%
  # Calculate the number of days between the onset of symptoms and the absence
  mutate(days_symptom_absent = as.numeric(difftime(date_absence,
                                                   date_symptom,
                                                   units = "days")))

# Matching of absences due to viral symptoms to positive saliva samples
# Range is set to four days
absence_viral <- link_absences_to_pathogens(absence_viral, range_days = 4) %>%
  mutate(pathogen = coalesce(pathogen, "No pathogen"),
         pathogen = factor(pathogen,
                           levels = c(pathgs, "No pathogen")))

# Absences due to other reasons than sickness
absence_non_viral <- abs_sus %>%
  filter(reason != "sickness") %>%
  dplyr::select(date_absence, date_return, student_id)

# Check whether absences due to other reasons than sickness are close to test results, 
# maybe no symptoms were reported
# Range here is set to two days
# absence_non_viral <- link_absences_to_pathogens(absence_non_viral, range_days = 2) %>%
#   filter(!is.na(pathogen))
# -> No absences due to other reasons than sickness were matched to positive test results

# Dates of Mondays during the study period, corresponding to weekends were students
# could not be tested
Mondays <- data.frame(
  date = seq.Date(as.Date("2024-01-23"), as.Date("2024-03-10"), by = "day")
) %>%
  mutate(
    weekday = weekdays(date),
    date = as.character(date)
  ) %>%
  filter(weekday %in% c("Monday")) %>%
  dplyr::select(date)

# All mandatory absences (weekends, school-free week)
abs_mandatory <- data.frame(
  date = seq.Date(as.Date("2024-01-20"), as.Date("2024-03-10"), 
                  by = "day")) %>%
  mutate(
    weekday = weekdays(date),
    student_id = Inf,
    date = as.character(date)
  ) %>%
  filter(weekday %in% c("Saturday", "Sunday")) %>%
  mutate(id = rep(1:(nrow(.) / 2), each = 2)) %>%
  dcast(id + student_id ~ weekday, value.var = "date") %>%
  mutate(across(c(Saturday, Sunday), as.Date)) %>%
  bind_rows(
    data.frame(
      student_id = Inf,
      Saturday = as.Date("2024-02-03"),
      Sunday = as.Date("2024-02-11")
    )
  ) %>%
  mutate(
    date_absence = Saturday,
    date_return = Sunday,
  ) %>%
  dplyr::select(student_id, date_absence, date_return) %>%
  filter(!is.na(date_absence)) %>%
  filter(row_number() != 3,
         row_number() != 4)

## Self-reported symptoms data #------------------------------------------------
symptoms <- absence_viral %>%
  dplyr::select(pathogen, symp_fever, symp_chills, symp_limb_pain, 
                symp_loss_of_taste, symp_loss_of_smell, symp_fatigue, symp_cough, 
                symp_cold, symp_diarrhea, symp_sore_throat, symp_headache, 
                symp_shortness_of_breath, symp_stomach_pains) %>%
  filter(pathogen != "No pathogen") %>%
  group_by(pathogen) %>%
  summarise(
    fever = mean(symp_fever, na.rm = T),
    chills = mean(symp_chills, na.rm = T),
    limb_pain = mean(symp_limb_pain, na.rm = T),
    loss_of_taste = mean(symp_loss_of_taste, na.rm = T),
    loss_of_smell = mean(symp_loss_of_smell, na.rm = T),
    fatigue = mean(symp_fatigue, na.rm = T),
    cough = mean(symp_cough, na.rm = T),
    cold = mean(symp_cold, na.rm = T),
    diarrhea = mean(symp_diarrhea, na.rm = T),
    sore_throat = mean(symp_sore_throat, na.rm = T),
    headache = mean(symp_headache, na.rm = T),
    shortness_of_breath = mean(symp_shortness_of_breath, na.rm = T),
    stomach_pains = mean(symp_stomach_pains, na.rm = T)
  )

## Social survey data #---------------------------------------------------------
pos_id <- unique(viral_loads$student_id)
soc_data <- soc_res %>%
  filter(time == "Baseline") %>%
  mutate(anyPos = factor(ifelse(student_id %in% pos_id, 1, 0), 
                         levels = c(0,1)),
         Club = factor(club_activity, levels = c("no", "yes"),
                       labels = c("Not active in a club", "Active in a club")),
         Gender =  factor(sex, levels = c("female", "male"),
                          labels = c("Female", "Male")),
         Siblings = factor(siblings, levels = c("no", "yes"),
                           labels = c("No siblings", "At least one sibling")),
         Class = factor(substr(student_id, 1,1),
                        levels = c("A", "B", "C", "D"),
                        labels = c("Class 1", "Class 2",
                                   "Class 3", "Class 4")))

## Construction of Positivity Periods #-----------------------------------------
# Include "time_id" of saliva test in "test_res"
test_res <- test_res %>%
  mutate(time_id = time_id[match(date, names(time_id))],
         # code if test was refused
         test_refused = is.na(pathogen_1) & absent == FALSE)

# Construction of positivity periods based on positive tests
infectious_periods <- viral_loads %>%
  dplyr::select(-panel, - class) %>%
  arrange(student_id, pathogen, date)

range.tests <- 1 # Number of tests between two positive tests allowed
ip <- 1 # Initial infectious period identifier
infectious_periods$infectious_period <- c(ip, rep(NA, nrow(infectious_periods)-1))
for (i in 2:nrow(infectious_periods)) {
  
  # Extract "time_id" where index student refused a test from "test_res"
  refused_time_id <- test_res %>%
    filter(student_id == infectious_periods$student_id[i],
           test_refused == TRUE) %>%
    pull(time_id)
  
  if (# Same pathogen
    infectious_periods$pathogen[i] == infectious_periods$pathogen[i - 1] &
    # Same student
    infectious_periods$student_id[i] == infectious_periods$student_id[i - 1] &
    # Time between positive tests is less or equal than a week
    difftime(infectious_periods$date[i], infectious_periods$date[i - 1],
             units = "days") <= 7 &
    # Not more than one test without positive test between two positive tests
    # Number of possible tests between two positive tests
    infectious_periods$time_id[i] - infectious_periods$time_id[i - 1] - 1 -
    # Number of tests refused during this period
    sum(refused_time_id > infectious_periods$time_id[i - 1] &
        refused_time_id < infectious_periods$time_id[i])
    # needs to be smaller than 1 test
    <= range.tests){
    # If condidtions fulfilled, 'i' and 'i-1' belong to the same infectious period
    infectious_periods$infectious_period[i] <- ip
    
    # Else, consider it a new infectious period
  } else {
    ip <- ip + 1
    infectious_periods$infectious_period[i] <- ip
  }
}

# Manual correction for student D5 for pathogen IFA
infectious_periods[infectious_periods$student_id == "D5" &
                     infectious_periods$date == "2024-01-31" &
                     infectious_periods$pathogen == "IFA", 
                   "infectious_period"] <- 78  

# Manual correction for student D5 for pathogen RSV
infectious_periods[infectious_periods$student_id == "D5" &
                     infectious_periods$date == "2024-01-31" &
                     infectious_periods$pathogen == "RSV", 
                   "infectious_period"] <- 83

# Manual correction for student D5 for pathogen AdV
infectious_periods[infectious_periods$student_id == "D5" &
                     infectious_periods$date == "2024-01-24" &
                     infectious_periods$pathogen == "AdV", 
                   "infectious_period"] <- 81

# Manual correction for D17
infectious_periods[infectious_periods$student_id == "D17" &
                     infectious_periods$pathogen == "IFB" &
                     infectious_periods$infectious_period == 73, 
                   "infectious_period"] <- 74

# Manual correction for D6
infectious_periods[infectious_periods$student_id == "D6" &
                     infectious_periods$pathogen == "RSV" &
                     infectious_periods$infectious_period == 86, 
                   "infectious_period"] <- 87


# Expansion of positivity periods based on absences
# abs_forced: mandatory absences (weekend, school-free week) for each student
absences_forced <- data.frame(
  student_id = rep(unique(sal_res$student_id), each = nrow(abs_mandatory)),
  date_absence = rep(abs_mandatory$date_absence,
                     length(unique(sal_res$student_id))),
  date_return = rep(abs_mandatory$date_return, 
                    length(unique(sal_res$student_id)))) %>%
  mutate(date_absence = as.Date(date_absence),
         date_return = as.Date(date_return))

# Check if mandatory absences are close to positive test results
# Conservative approach: 
# Range of one day (-> absence immediately before or after a positive test)
abs_forced <- link_absences_to_pathogens(absences_forced, 1) %>%
  mutate("date_symptom" = NA,
         type = "forced") %>%
  filter(!is.na(pathogen)) %>%
  arrange(student_id, pathogen, date_absence)

# "absences" contains all absences matched to positive test results (viral and
# forced absences)
absences <- bind_rows(absence_viral %>% 
                        mutate(date_return = as.Date(date_return), 
                               date_absence = as.Date(date_absence), 
                               type = "viral") %>%
                        filter(!is.na(test_id)) %>%
                        dplyr::select(student_id, date_absence, date_return, test_id,
                                      date_symptom, pathogen, type), abs_forced %>%
                        filter(!is.na(test_id)))

# Combine overlapping and sequential absence periods
# The following code merges absence periods assigned to the same
# positivity period, if they are overlapping or sequential
# Assign positivity period based on 'test_id' and 'pathogen' to absences
abs_possible <- absences %>%
  left_join(infectious_periods %>%
              dplyr::select(test_id, infectious_period, pathogen) %>%
              distinct(),
            by = c("test_id", "pathogen")) %>%
  mutate(absence_time = as.numeric(difftime(date_return, 
                                            date_absence, 
                                            units = "days"))) %>%
  # 1. Within each positivity period:
  #    a. If absence periods start on the same day, keep the longer one
  group_by(infectious_period, date_absence) %>%
  filter(absence_time == max(absence_time)) %>%
  ungroup() %>%
  #   b. If absence periods end on the same day, keep the longer one
  group_by(infectious_period, date_return) %>%
  filter(absence_time == max(absence_time)) %>%
  ungroup() %>%
  #   c. If absence periods are consecutive, merge them
  group_by(infectious_period) %>%
  arrange(date_absence) %>% 
  mutate(is_consecutive = lead(date_absence) == date_return + 1) %>%
  # If is_consecutive is TRUE, propagate the latest date_return over consecutive periods
  mutate(date_return = as.Date(cummax(
    as.numeric(if_else(is.na(is_consecutive) | !is_consecutive,
                       date_return, lead(date_return)))))) %>%
  filter(lag(is_consecutive) != TRUE | is.na(lag(is_consecutive))) %>%
  arrange(infectious_period) %>%
  #    c. If absence period x of a infectious period is contained in another absence 
  # period y of the same infectious period, remove x
  group_by(infectious_period) %>%
  arrange(date_absence) %>%
  mutate(is_contained = date_absence >= lag(date_absence) & 
           date_return <= lag(date_return)) %>%
  filter(is_contained == FALSE | is.na(is_contained)) %>%
  ungroup() %>%
  # Repeat using lead instead of lag
  group_by(infectious_period) %>%
  arrange(desc(date_absence)) %>%
  mutate(is_contained = date_absence <= lead(date_absence) & date_return >= lead(date_return)) %>%
  filter(is_contained == FALSE | is.na(is_contained)) %>%
  ungroup() %>%
  #     d. If absence periods are overlapping, merge them
  group_by(infectious_period) %>%
  arrange(date_absence) %>%
  mutate(is_overlapping = date_absence <= lag(date_return)) %>%
  # If is_overlapping is TRUE, propagate the latest date_return over overlapping periods
  mutate(date_return = as.Date(cummax(
    as.numeric(if_else(is.na(is_overlapping) | !is_overlapping,
                       date_return, lag(date_return)))))) %>%
  filter(lag(is_overlapping) != TRUE | is.na(lag(is_overlapping))) %>%
  arrange(infectious_period) %>%
  # Repeat using lead instead of lag
  group_by(infectious_period) %>%
  arrange(desc(date_absence)) %>%
  mutate(is_overlapping = date_absence <= lead(date_return)) %>%
  filter(is_overlapping == FALSE | is.na(is_overlapping)) %>%
  ungroup() %>%
  # Select relevant columns
  dplyr::select(-is_consecutive, -is_overlapping, -is_contained, -type) %>%
  # Recalculate absence time
  mutate(absence_time = as.numeric(difftime(date_return, 
                                            date_absence, 
                                            units = "days")),
         # Re-assign absence type (viral, forced)
         type_absence = if_else(date_return %in% abs_mandatory$date_return &
                                  date_absence %in% abs_mandatory$date_absence,
                                "forced", "viral"))

# Calculate "certain" and "uncertain" (censored) length of positivity periods

# Central objects: 
#   a. infectious_periods: contains the positivity periods based on positive tests
#   b. abs_possible: contains the absences matched to positivity periods

# Target variables:
#   a. observed: length of positivity period based on positive tests
#   b. possible: length of positivity period based on absences

# Target matrix
IC_time <- matrix(NA, nrow = length(unique(infectious_periods$infectious_period)),
                  ncol = 3)
IC_time[,1] <- unique(infectious_periods$infectious_period)
colnames(IC_time) <- c("infectious_period", "observed", "possible")

# Fill the target matrix
for (i in 1:nrow(IC_time)){
  # Infectious period
  ip <- IC_time[i,1]
  
  # Associated positive tests
  tests <- infectious_periods %>%
    filter(infectious_period == ip) %>%
    dplyr::select(date, time_id, test_id) %>%
    arrange(date)
  
  # Associated absences
  abs <- abs_possible %>%
    filter(infectious_period == ip) %>%
    dplyr::select(date_absence, date_return) %>%
    arrange(date_absence)
  
  # Check if one of the absences lays between one pair of positive tests
  if (nrow(tests) > 1 & nrow(abs) > 0){
    
    for (j in 1:(nrow(tests) -1)){
      # Remove absences that lay between two positive tests
      abs <- abs %>%
        filter(!(abs$date_absence >= tests$date[j] & 
                   abs$date_return <= tests$date[j + 1]))
    }
  }
  
  # Calculate observed length of positivity period
  IC_time[i,2] <- as.numeric(difftime(max(tests$date), 
                                      min(tests$date), units = "days")) + 1
  
  # Calculate the sum of all absence periods in days
  IC_time[i,3] <- sum(as.numeric(difftime(abs$date_return, 
                                          abs$date_absence, 
                                          units = "days")) + 1)
}

# Add observed and censored length of positivity periods to infectious_periods
infectious_periods_IC <- infectious_periods %>%
  left_join(data.frame(IC_time), by = "infectious_period") %>%
  mutate(interval_censored = possible > 0,
         interval_censored = if_else(is.na(interval_censored), 0, 1)) %>%
  # Add type of absence (viral, forced)
  left_join(abs_possible %>%
              dplyr::select(infectious_period, type_absence) %>%
              group_by(infectious_period) %>%
              arrange(desc(type_absence)) %>%
              filter(row_number() == 1),
            by = "infectious_period") %>%
  # Recode interval_censored
  mutate(interval_censored = if_else(possible > 0, 1, 0))

# Manual adjustment of type for B11, RSV
infectious_periods_IC[infectious_periods_IC$student_id == "B11" &
                        infectious_periods_IC$pathogen == "RSV" &
                        infectious_periods_IC$infectious_period == 17,
                      "type_absence"] <- "viral"

# Manual adjustment of type for B20, HRV 
infectious_periods_IC[infectious_periods_IC$student_id == "B20" &
                        infectious_periods_IC$pathogen == "HRV" &
                        infectious_periods_IC$infectious_period == 28,
                      "possible"] <- 0

# Manual adjustment of type for B3, HRV
infectious_periods_IC[infectious_periods_IC$student_id == "B3" &
                        infectious_periods_IC$pathogen == "HRV" &
                        infectious_periods_IC$infectious_period == 32,
                      "possible"] <- 9

# Manual adjustment of type for B9, HRV
infectious_periods_IC[infectious_periods_IC$student_id == "B9" &
                        infectious_periods_IC$pathogen == "HRV" &
                        infectious_periods_IC$infectious_period == 44,
                      "possible"] <- 9
##------------------------------------------------------------------------------

