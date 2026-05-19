# Outdoor Thermal Comfort and Behavioral Adaptation Data — Xixin Village

This repository contains the dataset supporting the manuscript:

> Gao P., Abd Ghafar A., Liu Z., Ghazalli A. J., Yeo L. B., & Wang L. 
> *Behavioral adaptation and outdoor thermal comfort of the rural elderly 
> in central China during the transitional season.* 

## Contents

| File | Description |
|------|-------------|
| `questionnaire_data.csv` | 196 valid responses with paired microclimate readings |
| `envimet_forcing_apr20.csv` | Hourly forcing data from Zhengzhou station (WMO 57083) for ENVI-met simulation |
| `data_dictionary.md` | Variable definitions, units, and coding schemes |
| `README.md` | This file |

## Data dictionary

Key variable groupings:

- **Microclimate**: Tg, Ta, RH, Va, Tmrt, PET
- **Demographics**: sex (1=male, 2=female), age (1=60–69, 2=70–79, 3=≥80)
- **Anthropometrics**: weight (1=<40kg, ..., 5=≥70kg), height (1=<150cm, ..., 5=≥180cm) — collected as categorical bands
- **Physiological**: activity met (W/m²), clo (clothing insulation)
- **Thermal evaluation**: TSV (−3 to +3), TCV (−2 to +2), TPV per variable, TAV (0/1)
- **Behavioral**: outdoor frequency, duration, location preferences, landscape preferences

## Field campaign

- **Site**: Xixin Village, Zhengzhou suburbs, Henan Province, China
- **Period**: 20–25 April 2025 (transitional season)
- **Sample**: 196 elderly residents (≥60 years), intercept survey at four outdoor sites (EW street, NS street, green land, fitness square)

## Citation

If you use this dataset, please cite the manuscript above.

## Contact

author: aygaopeng@gmail.com
