################################################################################
# Title: "data_preparation_descr_analysis.R"                                   #
# Description: Data preparation for the descriptive analysis of data           #
#              collected in the MCID_3 study.                                  #
################################################################################

## Dates of saliva sample collection
##------------------------------------------------------------------------------
test.dates <- unique(sal_res$date)[order(as.Date(unique(sal_res$date)))]
time_id <- 1:18
names(time_id) <- test.dates

## Saliva sample test results
##------------------------------------------------------------------------------
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

## Self-reported symptoms data
##------------------------------------------------------------------------------
symptoms <- unnest(infections_imp, cols = "symptoms") %>%
  dplyr::select(pathogen, symp_fever, symp_chills, symp_limb_pain, 
                symp_loss_of_taste, symp_loss_of_smell, symp_fatigue, symp_cough, 
                symp_cold, symp_diarrhea, symp_sore_throat, symp_headache, 
                symp_shortness_of_breath, symp_stomach_pains) %>%
  filter(!is.na(pathogen)) %>%
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

## Social survey data
##------------------------------------------------------------------------------
soc_data <- soc_res %>%
  filter(time == "Baseline") %>%
  mutate(Club = factor(club_activity, levels = c("no", "yes"),
                       labels = c("Not active in a club", "Active in a club")),
         Gender =  factor(sex, levels = c("female", "male"),
                          labels = c("Female", "Male")),
         Siblings = factor(siblings, levels = c("no", "yes"),
                           labels = c("No siblings", "At least one sibling")),
         Class = factor(substr(student_id, 1,1),
                        levels = c("A", "B", "C", "D"),
                        labels = c("Class 1", "Class 2",
                                   "Class 3", "Class 4")))


## School absenteeism data 
##------------------------------------------------------------------------------
# Absences due to sickness with respiratory symptoms
# absence_viral <- abs_sus %>%
#   filter(reason == "sickness", symp_respir == 1) %>%
#   # Calculate the number of days between the onset of symptoms and the absence
#   mutate(days_symptom_absent = as.numeric(difftime(date_absence,
#                                                    date_symptom,
#                                                    units = "days")))
# 
# # Matching of absences due to viral symptoms to positive saliva samples
# # Range is set to four days
# absence_viral <- link_absences_to_pathogens(absence_viral, range_days = 4) %>%
#   mutate(pathogen = coalesce(pathogen, "No pathogen"),
#          pathogen = factor(pathogen,
#                            levels = c(pathgs, "No pathogen")))

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

# Dates of Mondays/Fridays during the study period
Mondays <- data.frame(
  date = seq.Date(as.Date("2024-01-23"), as.Date("2024-03-10"), by = "day")
) %>%
  mutate(
    weekday = weekdays(date),
    date = as.character(date)
  ) %>%
  filter(weekday %in% c("Monday")) %>%
  dplyr::select(date)
Fridays <- as.Date(Mondays$date) -2

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

# Absences which are not due to viral symptoms (mandatory absences, other reasons)
absences_nonvirus <- data.frame(
  student_id = rep(unique(sal_res$student_id), each = nrow(abs_mandatory)),
  date_absent = rep(abs_mandatory$date_absence,
                     length(unique(sal_res$student_id))),
  date_return = rep(abs_mandatory$date_return, 
                    length(unique(sal_res$student_id)))) %>%
  mutate(date_absent = as.Date(date_absent),
         date_return = as.Date(date_return)) %>% 
  bind_rows(absence_non_viral %>%
              mutate(
                date_absence = as.Date(date_absence),
                date_return = as.Date(date_return)
              ))

## Inclusion of forced absences into absence data
infections_pos <- infections_imp %>%
  filter(!is.na(infection_id)) %>%
  dplyr::select(student_id, infection_id, pathogen, date_sample_first, date_sample_last, date_absent,
                date_return) %>%
  bind_rows(
    infections_imp %>% 
      filter(!is.na(infection_id)) %>%
      group_by(infection_id) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(student_id, infection_id, pathogen, date_sample_first, date_sample_last) %>%
      left_join(absences_nonvirus, by = "student_id", relationship = "many-to-many") %>%
      filter(
        date_sample_first == date_return + 1 | 
          date_sample_last == date_absence - 1
      )
  )

## Construction of positivity periods
##------------------------------------------------------------------------------
pos_periods <- infections_pos %>%
  mutate(
    absence_time = as.numeric(difftime(date_return, date_absent, units = "days"))
  ) %>%
  group_by(infection_id) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  # The following code is relevant for the construction of positivity periods when
  # multiple absences are matched to the same positivity period, hence these
  # are treated separately
  bind_rows(
    infections_pos %>%
    group_by(infection_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(
      absence_time = as.numeric(difftime(date_return, date_absent, units = "days")),
    ) %>%
    # 1. Within each positivity period:
    #    a. If absence periods start on the same day, keep the longer one
    group_by(infection_id, date_absent) %>%
    filter(absence_time == max(absence_time)) %>%
    ungroup() %>%
    #   b. If absence periods end on the same day, keep the longer one
    group_by(infection_id, date_return) %>%
    filter(absence_time == max(absence_time)) %>%
    ungroup() %>%
    #   c. If absence periods are consecutive, merge them
    group_by(infection_id) %>%
    arrange(date_absent) %>% 
    mutate(is_consecutive = lead(date_absent) == date_return + 1) %>%
    # If is_consecutive is TRUE, propagate the latest date_return over consecutive periods
    mutate(date_return = as.Date(cummax(
      as.numeric(if_else(is.na(is_consecutive) | !is_consecutive,
                         date_return, lead(date_return)))))) %>%
    filter(lag(is_consecutive) != TRUE | is.na(lag(is_consecutive))) %>%
    arrange(infection_id) %>%
    #    c. If absence period x of a infectious period is contained in another absence 
    # period y of the same infectious period, remove x
    group_by(infection_id) %>%
    arrange(date_absent) %>%
    mutate(is_contained = date_absent >= lag(date_absent) & 
             date_return <= lag(date_return)) %>%
    filter(is_contained == FALSE | is.na(is_contained)) %>%
    ungroup() %>%
    # Repeat using lead instead of lag
    group_by(infection_id) %>%
    arrange(desc(date_absent)) %>%
    mutate(is_contained = date_absent <= lead(date_absent) & date_return >= lead(date_return)) %>%
    filter(is_contained == FALSE | is.na(is_contained)) %>%
    ungroup() %>%
    #     d. If absence periods are overlapping, merge them
    group_by(infection_id) %>%
    arrange(date_absent) %>%
    mutate(is_overlapping = date_absent <= lag(date_return)) %>%
    # If is_overlapping is TRUE, propagate the latest date_return over overlapping periods
    mutate(date_return = as.Date(cummax(
      as.numeric(if_else(is.na(is_overlapping) | !is_overlapping,
                         date_return, lag(date_return)))))) %>%
    filter(lag(is_overlapping) != TRUE | is.na(lag(is_overlapping))) %>%
    arrange(infection_id) %>%
    # Repeat using lead instead of lag
    group_by(infection_id) %>%
    arrange(desc(date_absent)) %>%
    mutate(is_overlapping = date_absent <= lead(date_return)) %>%
    filter(is_overlapping == FALSE | is.na(is_overlapping)) %>%
    ungroup() %>%
    # Select relevant columns
    dplyr::select(-is_consecutive, -is_overlapping, -is_contained) %>%
    # Recalculate absence time
    mutate(absence_time = as.numeric(difftime(date_return, 
                                              date_absent, 
                                              units = "days"))
           )
  ) %>%
  # Calculate "certain" and "uncertain" (censored) length of positivity periods
  group_by(infection_id) %>%
  mutate(
    observed = difftime(date_sample_last, date_sample_first, units = "days") + 1,
    possible = case_when(
      is.na(absence_time) ~ 0,
      TRUE ~ absence_time
    ),
    possible = sum(absence_time) + 1
  ) %>%
  slice(1) %>%
  ungroup()



























