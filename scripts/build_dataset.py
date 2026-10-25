# -*- coding: utf-8 -*-
"""从原始提取文件重建分析数据集（COPD + Sepsis-3 + 入ICU后首次化验）

输出：
  cohort_72h.csv / cohort_24h.csv   每个 stay 一行
  列：人口学、生命体征、化验、评分、合并症、干预、结局
"""
import os
import pandas as pd, numpy as np, os, glob

B = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(B, "重建", "data")
os.makedirs(OUT, exist_ok=True)

def rd(name, cols=None):
    df = pd.read_csv(os.path.join(B, name), encoding="utf-8-sig", low_memory=False)
    df = df.drop_duplicates("stay_id")
    if cols:
        df = df[["stay_id"] + [c for c in cols if c in df.columns]]
    return df

base = pd.read_csv(os.path.join(B, "base_min_icu_mtepbq.csv"), encoding="utf-8-sig", low_memory=False)
d = base[["subject_id","hadm_id","stay_id","gender","dod","admittime","dischtime","los_hospital",
          "admission_age","race","hospital_expire_flag","icu_intime","icu_outtime","los_icu",
          "sofa_score","hosp_time_28d","hosp_outcome_28d"]].copy()

# ---- 合并各部分 ----
spec = {
 "copd":  ("metric_none_d_mtepbq_m2_1.csv",  ["chronic_obstructive_pulmonary_disease"]),
 "w":     ("metric_none_d_mtepbq_m3_1.csv",  ["weight_idx1"]),
 "h":     ("metric_none_d_mtepbq_m3_2.csv",  ["height_idx1"]),
 "sbp":   ("metric_none_d_mtepbq_m3_3.csv",  ["sbp_idx1"]),
 "dbp":   ("metric_none_d_mtepbq_m3_4.csv",  ["dbp_idx1"]),
 "mbp":   ("metric_none_d_mtepbq_m3_5.csv",  ["mbp_idx1"]),
 "hr":    ("metric_none_d_mtepbq_m3_6.csv",  ["heart_rate_idx1"]),
 "rr":    ("metric_none_d_mtepbq_m3_7.csv",  ["resp_rate_idx1"]),
 "temp":  ("metric_none_d_mtepbq_m3_8.csv",  ["temperature_idx1"]),
 "spo2":  ("metric_none_d_mtepbq_m3_9.csv",  ["spo2_idx1"]),
 "glu":   ("metric_none_d_mtepbq_m3_10.csv", ["glucose_idx1"]),
 "pt":    ("metric_none_d_mtepbq_m4_1.csv",  ["pt_idx1"]),
 "inr":   ("metric_none_d_mtepbq_m4_2.csv",  ["inr_idx1"]),
 "crp":   ("metric_none_d_mtepbq_m4_3.csv",  ["crp_idx1"]),
 "alb":   ("metric_none_d_mtepbq_m4_4.csv",  ["albumin_idx1","charttimediff_albumin_idx1"]),
 "glob":  ("metric_none_d_mtepbq_m4_5.csv",  ["globulin_idx1"]),
 "cl":    ("metric_none_d_mtepbq_m4_6.csv",  ["chloride_idx1"]),
 "ag":    ("metric_none_d_mtepbq_m4_7.csv",  ["aniongap_idx1"]),
 "bun":   ("metric_none_d_mtepbq_m4_8.csv",  ["bun_idx1"]),
 "na":    ("metric_none_d_mtepbq_m4_9.csv",  ["sodium_idx1"]),
 "k":     ("metric_none_d_mtepbq_m4_10.csv", ["potassium_idx1"]),
 "cr":    ("metric_none_d_mtepbq_m4_12.csv", ["creatinine_idx1"]),
 "alt":   ("metric_none_d_mtepbq_m4_13.csv", ["alt_idx1"]),
 "alp":   ("metric_none_d_mtepbq_m4_14.csv", ["alp_idx1"]),
 "ast":   ("metric_none_d_mtepbq_m4_15.csv", ["ast_idx1"]),
 "ldh":   ("metric_none_d_mtepbq_m4_16.csv", ["ld_ldh_idx1"]),
 "plt":   ("metric_none_d_mtepbq_m4_19.csv", ["platelet_idx1"]),
 "rdw":   ("metric_none_d_mtepbq_m4_21.csv", ["rdw_idx1"]),
 "wbc":   ("metric_none_d_mtepbq_m4_22.csv", ["wbc_idx1"]),
 "mono":  ("metric_none_d_mtepbq_m4_23.csv", ["monocytes_abs_idx1"]),
 "lym":   ("metric_none_d_mtepbq_m4_24.csv", ["lymphocytes_abs_idx1","charttimediff_lymphocytes_abs_idx1"]),
 "eos":   ("metric_none_d_mtepbq_m4_25.csv", ["eosinophils_abs_idx1"]),
 "baso":  ("metric_none_d_mtepbq_m4_26.csv", ["basophils_abs_idx1"]),
 "pco2":  ("metric_none_d_mtepbq_m4_27.csv", ["pco2_idx1"]),
 "po2":   ("metric_none_d_mtepbq_m4_28.csv", ["po2_idx1"]),
 "be":    ("metric_none_d_mtepbq_m4_29.csv", ["baseexcess_idx1"]),
 "ph":    ("metric_none_d_mtepbq_m4_30.csv", ["ph_idx1"]),
 "hb":    ("metric_none_d_mtepbq_m4_31.csv", ["hemoglobin_idx1"]),
 "mi":    ("metric_none_d_mtepbq_m5_1.csv",  ["myocardial_infarct"]),
 "chf":   ("metric_none_d_mtepbq_m5_2.csv",  ["congestive_heart_failure"]),
 "pvd":   ("metric_none_d_mtepbq_m5_3.csv",  ["peripheral_vascular_disease"]),
 "cvd":   ("metric_none_d_mtepbq_m5_4.csv",  ["cerebrovascular_disease"]),
 "dm":    ("metric_none_d_mtepbq_m5_5.csv",  ["diabetes"]),
 "renal": ("metric_none_d_mtepbq_m5_6.csv",  ["renal_disease"]),
 "cancer":("metric_none_d_mtepbq_m5_7.csv",  ["malignant_cancer"]),
 "liver": ("metric_none_d_mtepbq_m5_8.csv",  ["severe_liver_disease"]),
 "aps":   ("metric_none_d_mtepbq_m6_1.csv",  ["apsiii_min","apsiii_count"]),
 "sofa1": ("metric_none_d_mtepbq_m7_1.csv",  ["sofa_idx1"]),
 "saps":  ("metric_none_d_mtepbq_m8_1.csv",  ["sapsii_idx1"]),
 "vent":  ("metric_none_d_mtepbq_m9_1.csv",  ["ventilator_flag","ventilator_hours"]),
 "crrt":  ("metric_none_d_mtepbq_m10_1.csv", ["crrt_flag"]),
 "aki":   ("metric_none_d_mtepbq_m11_1.csv", ["aki_flag","aki_stage"]),
 "htn":   ("metric_none_d_mtepbq_m12_1.csv", ["hypertensive"]),
 "del":   ("metric_none_d_mtepbq_m1_1.csv",  ["delirium_flag"]),
}
for k, (fn, cols) in spec.items():
    d = d.merge(rd(fn, cols), on="stay_id", how="left")

# ---- 队列 ----
d = d[d.chronic_obstructive_pulmonary_disease == 1].copy()
d["los_icu"] = pd.to_numeric(d.los_icu, errors="coerce")
d = d[d.los_icu > 1]                                        # ICU 停留 >24h
d["d_alb"] = pd.to_numeric(d.charttimediff_albumin_idx1, errors="coerce")
d["d_lym"] = pd.to_numeric(d.charttimediff_lymphocytes_abs_idx1, errors="coerce")

num = ["albumin_idx1","lymphocytes_abs_idx1","rdw_idx1","aniongap_idx1","platelet_idx1","wbc_idx1",
       "monocytes_abs_idx1","eosinophils_abs_idx1","basophils_abs_idx1","bicarbonate_idx1" if "bicarbonate_idx1" in d else "bun_idx1"]
d["platelet_idx1"] = pd.to_numeric(d.get("platelet_idx1"), errors="coerce")

# ---- 数值清洗（明确写进 Methods）----
def clean(s, lo, hi):
    s = pd.to_numeric(s, errors="coerce")
    return s.where((s >= lo) & (s <= hi))

d["alb"]  = clean(d.albumin_idx1, 1.0, 6.0)
d["lym"]  = clean(d.lymphocytes_abs_idx1, 0.01, 15.0)
d["rdw"]  = clean(d.rdw_idx1, 9.0, 40.0)
d["ag"]   = clean(d.aniongap_idx1, 0.0, 60.0)
d["plt"]  = clean(d.platelet_idx1, 1, 2000)
d["wbc"]  = clean(d.wbc_idx1, 0.1, 200)
d["mono"] = clean(d.monocytes_abs_idx1, 0.0, 10)
d["eos"]  = clean(d.eosinophils_abs_idx1, 0.0, 10)
d["baso"] = clean(d.basophils_abs_idx1, 0.0, 10)
d["neut"] = d.wbc - d[["lym","mono","eos","baso"]].sum(axis=1, min_count=1)
d.loc[(d.neut <= 0) | (d.neut > 100), "neut"] = np.nan

# ---- 指标 ----
d["PNI"]  = 10*d.alb + 5*d.lym
d["RAR"]  = d.rdw/d.alb
d["ACAG"] = d.ag + 0.25*(4.0 - d.alb)
d["NLR"]  = d.neut/d.lym
d["PLR"]  = d.plt/d.lym
d["SII"]  = d.plt*d.neut/d.lym
d["SIRI"] = d.neut*d.mono/d.lym
d["BMI"]  = pd.to_numeric(d.weight_idx1, errors="coerce")/(pd.to_numeric(d.height_idx1, errors="coerce")/100)**2

# ---- 结局 ----
for c in ["admittime","dischtime","dod","icu_intime"]:
    d[c] = pd.to_datetime(d[c], errors="coerce")
d["d_adm"] = (d.dod - d.admittime).dt.total_seconds()/86400
d["d_icu"] = (d.dod - d.icu_intime).dt.total_seconds()/86400    # 以入 ICU 为起点（主分析）
d["death90"] = ((d.d_icu >= 0) & (d.d_icu <= 90)).astype(int)
d["t90"] = np.where(d.death90 == 1, d.d_icu.clip(lower=0.04, upper=90), 90)
d["death90_adm"] = ((d.d_adm >= 0) & (d.d_adm <= 90)).astype(int)
d["death28"] = pd.to_numeric(d.hosp_outcome_28d, errors="coerce").fillna(0).astype(int)

d["male"] = (d.gender.astype(str).str.upper().str.startswith("M")).astype(int)
d["race_white"] = d.race.astype(str).str.upper().str.contains("WHITE").astype(int)

for w in (24, 48, 72):
    d[f"win{w}"] = ((d.d_alb >= 0) & (d.d_alb <= w) & (d.d_lym >= 0) & (d.d_lym <= w)).astype(int)

keep = ["stay_id","subject_id","hadm_id","win24","win48","win72","d_alb","d_lym",
        "PNI","RAR","ACAG","NLR","PLR","SII","SIRI","alb","lym","rdw","neut","plt","wbc","mono",
        "crp_idx1","hb","globulin_idx1","bun_idx1","creatinine_idx1","sodium_idx1","potassium_idx1",
        "chloride_idx1","bicarbonate_idx1","ph_idx1","pco2_idx1","po2_idx1","baseexcess_idx1",
        "ld_ldh_idx1","inr_idx1","pt_idx1","glucose_idx1",
        "sbp_idx1","dbp_idx1","mbp_idx1","heart_rate_idx1","resp_rate_idx1","temperature_idx1","spo2_idx1",
        "admission_age","male","race_white","BMI","los_hospital","los_icu",
        "sofa_score","sofa_idx1","sapsii_idx1","apsiii_min","apsiii_count",
        "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease","cerebrovascular_disease",
        "diabetes","renal_disease","malignant_cancer","severe_liver_disease","hypertensive",
        "ventilator_flag","ventilator_hours","crrt_flag","aki_flag","aki_stage","delirium_flag",
        "hospital_expire_flag","death28","death90","death90_adm","t90","hosp_time_28d","hosp_outcome_28d"]
keep = [c for c in keep if c in d.columns]
d = d[keep]
d.to_csv(os.path.join(OUT, "cohort_all.csv"), index=False, encoding="utf-8-sig")

print("COPD+Sepsis3+ICU>24h:", len(d))
for w in (72, 48, 24):
    s = d[d[f"win{w}"] == 1]
    print(f"  win{w}: n={len(s)}, 90d死亡={int(s.death90.sum())}, 院内死亡={int(s.hospital_expire_flag.sum())}")
print("\n关键变量缺失率（72h 子集，n=%d）:" % int((d.win72==1).sum()))
s = d[d.win72 == 1]
for c in ["PNI","RAR","ACAG","NLR","SII","SIRI","sapsii_idx1","sofa_score","BMI","crp_idx1","rdw","alb","lym"]:
    print(f"  {c:18s} {100*s[c].isna().mean():5.1f}%")
