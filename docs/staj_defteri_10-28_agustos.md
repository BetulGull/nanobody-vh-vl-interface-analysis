# Staj Defteri — 10 Ağustos – 28 Ağustos 2026

**Proje:** VHH (tek domainli deve antikoru) yapısal mekanizma analiz hattı

**Teknik özet:** SAbDab kaynaklı ~4.500 deneysel protein yapısı (mmCIF) üzerinde çalışan,
aşamalı (Stage 0–12), deterministik ve tekrar üretilebilir bir Python veri işleme +
istatistik + görselleştirme hattının geliştirilmesi.

**Projenin cevap aradığı soru:** Klasik antikorlarda ağır zincir değişken domaini (VH),
hafif zincir domaini (VL) ile eşleşir ve VL, VH yüzeyinin bir bölümünü fiziksel olarak
kapatır. VHH domainleri ise VL olmadan tek başına çalışır. Soru: VL'nin kapattığı bu yüzey
("former interface" / eski arayüz) VHH'de nasıl farklılaşmıştır? Bu, yazılım tarafında şu
probleme indirgenmiştir: **iki kohortu (VHH vs klasik VH) aynı koordinat sisteminde (IMGT
numaralandırması), aynı eşiklerle, aynı istatistiksel testlerle, hiçbir manuel seçim
yapmadan karşılaştıran bir pipeline yazmak.**

---

## 10 Ağustos, Pazartesi — Proje kurulumu, veri alımı ve dizin mimarisi

**Yapılan iş.** İlk gün proje iskeletini ve veri alım katmanını kurdum. Hesaplama
araçlarının bağımlılıkları çakıştığı için tek bir sanal ortam yerine üç izole `micromamba`
ortamı tanımladım (`environment-*.yml`): ana bilimsel ortam (Python 3.11), etkileşim analizi
ortamı (Python 3.9, eski bağımlılıklar) ve elektrostatik ortamı. Dış araçlar kütüphane
olarak import edilmiyor, `subprocess` ile çağrılıyor; böylece sürüm çakışması pipeline'ı
kırmıyor.

Dizin mimarisi katmanlı: `data/raw_extracted` (salt-okunur ham veri), `data/interim`,
`data/processed` (Parquet), `src/` (modüller), `configs/config.yaml` (tüm eşikler tek
yerde), `results/`, `logs/`, `tests/`. Kural: **hiçbir aşama ham veriyi değiştiremez**, her
aşama kendi çıktısını yazar. Böylece herhangi bir aşama baştan çalıştırılabiliyor.

Arşivlerin açılması için `safe_extract` yazdım. Üç problem çözüldü: (1) tar path-traversal
güvenlik açığı, (2) idempotency — aynı arşiv tekrar açılırsa var olan dosya SHA256 ile
karşılaştırılıyor, farklıysa üzerine yazılmıyor, log'lanıyor, (3) macOS APFS büyük/küçük
harf duyarsız olduğu için sadece zincir harfi farkıyla ayrışan dosyaların (`_B_` vs `_b_`)
çakışması; bunlar yeniden adlandırılıp bir manifest'e kaydediliyor.

```python
# src/ingestion/safe_extract.py
def _is_safe_member(member: tarfile.TarInfo, dest_root: Path) -> bool:
    if member.name.startswith("/"):
        return False
    dest_path = (dest_root / member.name).resolve()
    try:
        dest_path.relative_to(dest_root.resolve())    # path-traversal koruması
    except ValueError:
        return False
    return True

def _actual_case_on_disk(path: Path):
    """APFS/HFS+ case-insensitive olduğu için sadece harf büyüklüğüyle
    ayrışan arşiv üyeleri çakışır; gerçek disk adını döndürür."""
    target_lower = path.name.lower()
    for entry in os.listdir(path.parent):
        if entry.lower() == target_lower:
            return entry
    return None
```

**Çıktı.** Figür değil: `extraction_manifest.tsv`, SHA256 provenance tablosu ve kurulum
log'ları (`logs/setup_commands.log`, `logs/extraction_conflicts.log`).

**Yorum.** Bu gün bilimsel bir sonuç üretmedi; tüm sonraki aşamaların üzerine kurulacağı
değişmez (immutable) veri katmanı ve izlenebilirlik altyapısı sağlandı.

---

## 11 Ağustos, Salı — Stage 1: metadata birleştirme ve kural tabanlı kohort sınıflandırma

**Yapılan iş.** Yapıların hangi kohorta ait olduğunu belirleyen sınıflandırma katmanını
yazdım. Girdi: birden fazla CSV export + arşiv manifesti. Bunları tek bir "master" tabloda
birleştirip alan normalizasyonu yaptım (`resolution` alanı serbest metin geliyordu, sayısala
zorlandı; deney yöntemi normalize edildi).

Sınıflandırma bilinçli olarak **kural tabanlı ve öncelik sıralı** yazıldı, ML kullanılmadı:
önce yapılandırılmış kanıt (`heavy_taxid` → tür kimliği), o yoksa tür adı metni, o da yoksa
serbest metin alanlarında anahtar kelime araması (mühendislik/hümanizasyon/sentetik/mutasyon
terimleri). Kritik nokta: her satır için sadece bir etiket değil, **karar gerekçesi de
string olarak** yazılıyor (`species_evidence_field`, `exclusion_reason`). Böylece 4.500
satırlık bir sınıflandırma manuel olarak denetlenebiliyor — bu bir denetim izi (audit trail)
tasarımı.

QC bayrakları da bu aşamada üretildi (çözünürlük ≤ 3 Å, yöntem, zincir bütünlüğü); dışlanan
yapılar silinmiyor, `False` bayrağıyla tabloda kalıyor ki huni (funnel) sayıları
raporlanabilsin.

```python
# src/qc/stage1_pipeline.py
CAMELID_TAXIDS = {...}; HUMAN_TAXID = 9606; SYNTHETIC_TAXID = 32630

def classify_sdh_row(row) -> dict:
    # Öncelik 1: yapılandırılmış taksonomi kimliği
    if heavy_taxid is not None:
        species_evidence_field = "heavy_taxid"
        camelid_species_flag = heavy_taxid in CAMELID_TAXIDS
    # Öncelik 2: tür adı metni  →  Öncelik 3: serbest metinde anahtar kelime
    eng_hits, _ = keyword_hits(text_fields, ENGINEERING_KEYWORDS)
    syn_hits, _ = keyword_hits(text_fields, SYNTHETIC_KEYWORDS)
    ...

df["qc_resolution_numeric"] = df["resolution"].apply(coerce_resolution)
df["qc_resolution_le_3"]    = df["qc_resolution_numeric"].apply(
    lambda r: r is not None and r <= 3.0)
```

**Çıktı.** Figür değil: Stage-1 çözümlenmiş metadata tablosu — **1.458 VHH + 3.092 eşleşmiş
(VH–VL) yapı**, her satırda gerekçe ve QC bayrakları.

**Yorum.** Kohortlar bu noktada ilk kez tanımlandı; sonraki tüm karşılaştırmaların örneklem
tanımı bu tablodan geliyor.

---

## 12 Ağustos, Çarşamba — Stage 2: yapı QC, temas arama algoritması ve birim testler

**Yapılan iş.** İki iş yaptım: mmCIF dosyalarının yapısal kalite kontrolü ve projenin en sık
kullanılan çekirdek algoritması olan **temas (contact) aramanın** yazılması.

Temas sorusu şu: iki atom kümesi arasındaki en küçük ağır-atom mesafesi bir eşiğin altında
mı? Naif çözüm her atom çiftini gezmek — O(N·M). Bunu referans (doğruluk kaynağı) olarak
yazdım, üretimde ise `scipy.spatial.cKDTree` ile uzamsal indeksleme kullandım; sorgu
maliyeti pratikte O(N·log M)'e düşüyor. İki algoritmanın **birebir aynı sonucu** vermesi bir
testle garanti altına alındı — optimize sürümün doğruluğu bu eşdeğerlikle kanıtlanıyor.

Ayrıca: alternatif konformasyon (altloc) politikası (en yüksek doluluk, eşitlikte
alfabetik), model seçimi, eksik atom tespiti ve her 100 yapıda bir checkpoint yazımı eklendi
— 4.500 yapılık bir çalıştırma yarıda kesilirse baştan başlanmıyor.

```python
# src/qc/contact_search.py
def optimized_contacts(self, cdr_atoms, antigen_atoms):
    if len(cdr_atoms) == 0 or len(antigen_atoms) == 0:
        return False, [], []
    tree = cKDTree(antigen_atoms)                             # O(M log M) inşa
    hits_a, hits_b = [], set()
    for i, cdr_coord in enumerate(cdr_atoms):
        nearby = tree.query_ball_point(cdr_coord, self.threshold)   # yarıçap sorgusu
        if nearby:
            hits_a.append(i); hits_b.update(nearby)
    return len(hits_a) > 0, sorted(hits_a), sorted(hits_b)

# tests/test_contact_search.py — eşik sınırı, boş dizi, brute-force eşdeğerliği
def test_optimized_matches_brute_force_large():
    assert optimized_contacts(A, B) == brute_force_contacts(A, B)
```

**Çıktı.** Figür değil: 10/10 geçen birim test (`10 passed in 0.32s`) ve yapı-QC bayrak
tablosu.

**Yorum.** Eşik değeri (4,5 Å) tek bir config dosyasından geliyor; sonraki tüm aşamalar aynı
temas tanımını çağırıyor, yani sonuçlar arasında tanım tutarsızlığı imkânsız.

---

## 13 Ağustos, Perşembe — Figure 1: eski arayüzün ampirik olarak tanımlanması

**Yapılan iş.** Projenin merkezî tanımını ürettim: "VL'nin kapattığı yüzey" varsayılmadı,
**ölçülerek** bulundu. Yöntem, her eşleşmiş yapı için erişilebilir yüzey alanının (SASA) iki
kez hesaplanması: bir kez sadece VH atomlarıyla, bir kez VH+VL birlikte. Farkı pozitif olan
kalıntılar VL tarafından kapatılıyor demektir.

```python
# src/qc/stage3_interface_derivation.py
per_res_alone,   _ = calc_sasa(vh_atoms)               # FreeSASA, Lee-Richards
per_res_with_vl, _ = calc_sasa(vh_atoms + vl_atoms)
delta = sasa_alone - sasa_with_vl                      # VL nedeniyle kapanan alan
interface_1  = delta >= 1.0
interface_5  = delta >= 5.0        # birincil eşik (Å²)
interface_10 = delta >= 10.0       # üçü birden yazılır → duyarlılık analizi bedava
```

Bu residü düzeyindeki sonuç, IMGT pozisyonlarına göre kohort genelinde toplandı. Bir
pozisyonun "arayüz pozisyonu" sayılması için kohortun en az %25'inde eşiği geçmesi gerekiyor:

```python
# src/qc/stage3_position_aggregation.py
PRIMARY_FREQUENCY_THRESHOLD = 0.25
df['primary_binary_interface'] = (
    (df['interface_frequency_5']  >= PRIMARY_FREQUENCY_THRESHOLD) &
    (df['observation_fraction']   >= PRIMARY_OBSERVATION_FRACTION_THRESHOLD)
)
```

Sonuç: **28 IMGT pozisyonluk ampirik eski arayüz kümesi**. Ayrıca aynı diziyi paylaşan
yapıların sonucu şişirmemesi için redundans kümeleme yapıldı ve her kümeden deterministik
bir temsilci seçildi (kapsam, sonra çözünürlük, sonra kimlik sırasına göre) → analiz kohortu
**533 VHH / 1.429 VH**.

**Çıktı.** `Figure1.png` — (A/B) arayüzün yapı üzerinde gösterimi, (C) kohort boyutları,
(D/E) IMGT hizalı dizi logoları.

**Yorum.** Arayüz artık literatürden alınan bir varsayım değil, veriden türetilmiş ve
dondurulmuş (frozen) bir pozisyon listesi; sonraki 10 figürün tamamı bu listeyi girdi olarak
kullanıyor.

---

## 14 Ağustos, Cuma — Figure 2: pozisyon bazlı dizi karşılaştırması ve ıraksama

**Yapılan iş.** 28 pozisyonun her birinde iki kohortun amino asit dağılımını çıkardım ve
farkı ölçtüm. Her pozisyon için 20 elemanlı bir olasılık vektörü üretiliyor; iki dağılımın
farkı **Jensen–Shannon ıraksaması (JSD)** ile ölçülüyor. JSD seçildi çünkü simetrik, sınırlı
(0–1) ve sıfır frekanslarda tanımsız olmuyor — KL ıraksamasının aksine. `scipy` mesafe
döndürdüğü için karesi alınıyor.

```python
# src/qc/stage4_position_comparison.py
counts = grp['amino_acid'].value_counts()
freqs  = np.array([counts.get(aa, 0) / n for aa in CANONICAL_20])   # 20'lik dağılım
...
jsd = float(jensenshannon(vhh_freqs, vh_freqs, base=2) ** 2)   # mesafe² = ıraksama
comparison['delta_mean_KD'] = comparison['mean_KD_VHH'] - comparison['mean_KD_VH']
```

Ayrıca pozisyon × amino asit boyutunda işaretli bir zenginleşme matrisi üretildi (VHH
frekansı − VH frekansı), böylece "fark var" ifadesi "hangi harf hangi yöne" düzeyine
indirgendi. Görselleştirmede logo panelleri, ısı haritası ve lolipop grafiği aynı x-ekseni
(IMGT pozisyonu) üzerinde hizalandı — okuyucunun üç paneli dikey olarak takip edebilmesi
için.

**Çıktı.** `Figure2.png` — (A/B) VHH ve VH dizi logoları, (C) kalıntı zenginleşme matrisi,
(D) pozisyon bazlı JSD.

**Yorum.** Fark tüm yüzeye yayılmıyor; birkaç pozisyonda (özellikle FR2 bölgesindeki
42/49/50/52 dörtlüsü) JSD ~0,7–0,8'e çıkarken diğer pozisyonlar neredeyse aynı — yani sinyal
yerel ve yoğunlaşmış.

---

## 17 Ağustos, Pazartesi — Figure 3: CDR3 uzunluğu × grup etkileşim modeli

**Yapılan iş.** Bu güne kadar üretilen ölçümler yapı düzeyinde bir tabloda toplanmıştı;
artık "iki grup arasında fark var mı" sorusunu **istatistiksel olarak** sormak gerekiyordu.
Yaygın hata, her grup için ayrı regresyon çalıştırıp iki p-değerini karşılaştırmak. Bunun
yerine tek bir önceden ilan edilmiş (predeclared) model kuruldu; grup farkı doğrudan
etkileşim katsayısı β₃ ile test ediliyor:

Y = β₀ + β₁·CDR3_uzunluk + β₂·grup + β₃·(CDR3_uzunluk × grup)

```python
# src/qc/final_cdr3_length_group_interactions.py
def fit_interaction(d, outcome):
    sub = d[[outcome, 'CDR3_length', 'group_bin']].dropna().copy()
    sub['length_x_group'] = sub['CDR3_length'] * sub['group_bin']   # group_bin: VHH=1, VH=0
    X = sm.add_constant(sub[['CDR3_length', 'group_bin', 'length_x_group']])
    model = sm.OLS(sub[outcome], X).fit(cov_type='HC3')   # heteroskedastisiteye dayanıklı SE
    return {'slope_VH':  model.params['CDR3_length'],
            'slope_VHH': model.params['CDR3_length'] + model.params['length_x_group'],
            'beta_interaction': model.params['length_x_group'],
            'CI_low': model.conf_int().loc['length_x_group'][0],
            'p_interaction': model.pvalues['length_x_group']}
```

Varyans uzunlukla arttığı için HC3 dayanıklı standart hatalar kullanıldı. Çoklu test problemi
için sonuçlar **altı önceden tanımlı aileye** bölünüp her aile içinde ayrı ayrı BH-FDR
düzeltmesi yapıldı (aileler asla birleştirilmedi) — bu, testleri sonradan seçip anlamlı
çıkanı raporlama riskini kapatıyor.

**Çıktı.** `Figure3.png` — üç örnek sonuç: CDR3 yükü, arayüz hidropati ve çerçeve teması;
alt panelde uzunluk gruplarına göre pozisyon bazlı hidropati ısı haritası.

**Yorum.** İki grubun regresyon doğrularının eğimleri belirgin biçimde farklı (β₃ anlamlı,
q < 0,001); yani ilişki VHH'ye özgü, genel bir antikor eğilimi değil.

---

## 18 Ağustos, Salı — **Raporlu (izinli)**

---

## 19 Ağustos, Çarşamba — Figure 4A: CDR3 ile eski arayüz arasındaki temas geometrisi

**Yapılan iş.** Stage-2'de yazdığım temas mantığını yeni bir soruya uyguladım: domainin
kendi CDR3 ilmeği, eski arayüzün 15 pozisyonluk alt kümesine ne kadar temas ediyor? Pipeline
önce **uzun (long) formatta** ham mesafe tablosu üretiyor — her (yapı, CDR3 pozisyonu, hedef
pozisyon) üçlüsü için tek satır; eşikleme ve olasılığa çevirme sonra, ayrı bir adımda
yapılıyor. Böylece 4,0 / 4,5 / 5,0 Å duyarlılık analizi tekrar tarama gerektirmiyor.

```python
# src/qc/figure4_geometry_pipeline.py
def min_dist(atoms1, atoms2):                 # vektörize ağır-atom minimum mesafesi
    p1 = np.array([[a['x'], a['y'], a['z']] for a in atoms1])
    p2 = np.array([[a['x'], a['y'], a['z']] for a in atoms2])
    diff = p1[:, None, :] - p2[None, :, :]
    return float(np.sqrt((diff ** 2).sum(axis=-1)).min())

rows.append({'INSTANCE': instance, 'group': group,
             'cdr3_position': cdr3_pos, 'target_position': tpos,
             'min_heavy_atom_distance': d, 'tetrad': tpos in TETRAD})
```

Kritik yazılımsal detay paydada: bazı yapılarda bazı pozisyonlar hiç gözlenmiyor. Bunları 0
kabul etmek olasılığı yapay olarak düşürürdü; gözlenmeyen çift matrise hiç girmiyor ve NaN
kalıyor. Bu kural bir testle kilitlendi:

```python
# tests/test_figure4A.py
def test_missing_pairs_excluded_not_zero_filled():
    denom = pd.read_csv(SRC / "Figure4A_contact_denominators.csv")
    p     = pd.read_csv(SRC / "Figure4A_VH_contact_probability.csv", index_col=0)
    zero_denom = denom[(denom.group == 'VH') & (denom.N_evaluable == 0)]
    for _, row in zero_denom.iterrows():
        assert pd.isna(p.loc[row['target_position'], row['cdr3_label']])
```

Grup farkı için dağılımlar normal olmadığından parametrik olmayan Mann-Whitney U + BH-FDR
kullanıldı. Yapı görselleri PyMOL ile betiklendi; **hiçbir örnek yapı gözle seçilmedi**, her
biri "grubun medyanına en yakın yapı" kuralıyla belirlendi ve seçim CSV'ye yazıldı:

```python
# scripts/figure4A_render_structural_context.py
cmd.color('iface_orange', f'{obj} and chain {sel_h} and resi {"+".join(map(str, SHIELDABLE_15))}')
cmd.color('cdr3_red',     f'{obj} and chain {sel_h} and resi {"+".join(map(str, CDR3_RANGE))}')
cmd.show('spheres',       f'{obj} and chain {sel_h} and resi {"+".join(map(str, TETRAD))} and name CA')
```

**Çıktı.** `Figure4A.png` — (A) yapısal bağlam, (B/C) grup bazlı temas olasılığı, (D)
işaretli fark matrisi, (E) yapı başına geometri özeti + anlamlılık.

**Yorum.** VHH'de CDR3, arayüzün daha geniş bir bölümüyle temas ediyor (ortalama 0,341 vs
0,267) ve arayüze daha yakın duruyor (7,95 Å vs 8,66 Å); fark matrisi bu farkın birkaç
pozisyonda toplandığını gösteriyor.

---

## 20 Ağustos, Perşembe — Figure 4B: literatürdeki bir analiz mantığının uyarlanması

**Yapılan iş.** Referans bir makalenin (Zhang 2025) analiz mantığını kendi
karşılaştırmamıza uyarladım. Buradaki yazılımsal disiplin şuydu: **referansın kendi eşiği
korundu** (4 Å), bizim projenin 4,5 Å eşiğiyle karıştırılmadı; iki eşik ayrı ayrı belgelendi.

Üç yeni ölçüm modülü yazıldı: (1) CDR3'ün diğer bölgelerle etkileşim sayımı, (2) domain SASA
ve yüzey hidrofobisitesi (SASA × Kyte-Doolittle), (3) DSSP ile ikincil yapı ataması. DSSP
harici bir ikili (binary); conda ortam değişkenleri altında segfault verdiği için
temizlenmiş bir ortamla `subprocess` üzerinden çağrıldı — dış araç entegrasyonunda tipik bir
sorun:

```python
# src/qc/figure4B_dssp.py
# mkdssp (homebrew binary) conda/micromamba ortam değişkenleriyle çağrıldığında segfault veriyor
_MKDSSP_BIN  = '/opt/homebrew/bin/mkdssp'
ZHANG_STATES = {'G', 'T', 'E', 'S', 'C'}     # 8 durumdan rapor alt kümesine indirgeme

result = subprocess.run([_MKDSSP_BIN, '--output-format', 'dssp', str(pdb_path), str(dssp_out)],
                        capture_output=True, text=True, env=_clean_env())
if result.returncode != 0:
    return [], f"mkdssp failed: {result.stderr[:200]}"    # sessiz başarısızlık yok
zhang_state = 'C' if ss == '-' else (ss if ss in ZHANG_STATES else 'OTHER')
```

Ayrıca iki duyarlılık analizi eklendi: ortak destek (common support — VH kohortu VHH'nin
uzunluk aralığına kısıtlanarak) ve uzunluğa normalize edilmiş metrik. Amaç, "uzun CDR3 daha
çok kalıntı içerdiği için daha çok temas eder" gibi trivial bir açıklamayı elemek.

**Çıktı.** `Figure4B.png` — (A) kısa/uzun CDR3 örnek yapıları, (B–D) CDR3'ün CDR1/FR2/CDR2
ile etkileşim olasılığı, (E) korelasyon + etkileşim testi özet tablosu.

**Yorum.** CDR3–FR2 etkileşiminin uzunlukla artışı belirgin şekilde VHH'ye özgü (r = 0,60 vs
0,19) ve normalizasyondan sonra da ayakta kalıyor; buna karşılık CDR3–CDR1 metriğinde
gruplar arası anlamlı fark **yok** — negatif sonuç da raporlandı.

---

## 21 Ağustos, Cuma — Figure 5: eski arayüzün hidrofobik yüzey skoru (PSH)

**Yapılan iş.** Yüzeyin kimyasal karakterini sayısallaştıran patch skorlama modülünü yazdım.
Formül literatürden tahmin edilmedi; referans aracın kaynak kodundan birebir çıkarıldı ve
her sabit yorum satırında "source-confirmed" olarak işaretlendi. Skor, yüzeyde birbirine
yakın kalıntı çiftleri üzerinden hesaplanan mesafe-ağırlıklı bir toplam:

```python
# src/qc/stage8_patch_scoring.py
EXPOSURE_THRESHOLD_PCT = 7.5   # source-confirmed: yan zincir bağıl erişilebilirlik >= %7,5
PATCH_PAIR_CUTOFF_A    = 7.5   # source-confirmed: çift mesafe kesimi

def compute_patch_scores(residues, vicinity):
    psh = ppc = pnc = 0.0
    for k1 in vicinity:
        for k2 in vicinity:                       # sıralı çiftler
            if k1 == k2: continue
            d = _min_heavy_atom_distance(residues[k1]['atoms'], residues[k2]['atoms'])
            if d > PATCH_PAIR_CUTOFF_A or d == 0: continue
            r1, r2 = residues[k1], residues[k2]
            psh += (r1['kyte_doolittle_norm'] * r2['kyte_doolittle_norm']) / (d * d)
            c1, c2 = r1['charge_class'], r2['charge_class']
            if   c1 > 0 and c2 > 0: ppc += (c1 * c2) / (d * d)
            elif c1 < 0 and c2 < 0: pnc += (c1 * c2) / (d * d)
    return {'PSH': psh, 'PPC': ppc, 'PNC': pnc}
```

Ayrıca bağlı bileşen (connected component) analizi: yüzeydeki hidrofobik kalıntılar 4 Å
komşuluk kriteriyle bir grafa dönüştürülüp kümelendi; en büyük bileşen "en büyük hidrofobik
yama" olarak raporlandı. Panel görselleri yine deterministik kuralla seçildi: grubun
medyanına ve 95. yüzdeliğine en yakın yapılar.

**Çıktı.** `Figure5.png` — (A) PSH dağılımları, (C/D) 95. yüzdelik yapıların yüzey
renklendirmesi, (E/H) en büyük bağlı hidrofobik yama (siyah), (F/G) medyan yapılar.

**Yorum.** VHH dağılımı daha düşük PSH'ye kaymış (medyan 76,0 vs 84,7); hem medyan hem 95.
yüzdelik yapılarda VHH'nin sürekli hidrofobik yaması VH'ninkinden küçük görünüyor.

---

## 24 Ağustos, Pazartesi — Figure 6: CDR3–arayüz paketlenmesi (gömülü yüzey alanı)

**Yapılan iş.** Temas *sayısı* yerine sürekli bir geometrik büyüklük olan **karşılıklı gömülü
yüzey alanını (BSA)** hesapladım. Tanım üç SASA çağrısının cebirsel farkı:

```python
# src/qc/stage7_bsa.py
bsa_cdr3_side   = sasa_a_alone - sum_a_in_complex
bsa_target_side = sasa_b_alone - sum_b_in_complex
bsa_pair_total  = sasa_a_alone + sasa_b_alone - sasa_ab_total

# tutarlılık kapısı: iki ayrı yoldan hesaplanan toplam eşleşmezse kalıntı üyeliği hatalıdır
consistency_ok = abs((bsa_cdr3_side + bsa_target_side) - bsa_pair_total) < 1e-6
return {'BSA_CDR3_side': bsa_cdr3_side, 'BSA_target_side': bsa_target_side,
        'BSA_pair_total': bsa_pair_total, 'bsa_consistency_ok': consistency_ok}
```

Buradaki kontrol bilimsel değil **yazılımsal** bir güvenlik ağı ve tüm kohort için
çalıştırıldı. Sonrasında 17 Ağustos'taki model formu bu çıktıya uygulandı (BSA ~ uzunluk +
grup + uzunluk×grup, HC3). Trivial açıklamayı elemek için iki bağımsız normalizasyon kontrolü
daha eklendi: kalıntı başına BSA ve sınırlı (0–1) paketlenme oranı.

**Çıktı.** `Figure6.png` — (A) CDR3 uzunluğu vs BSA regresyonu ve %95 güven bantları, (B/C)
Q25/Q75 uzunluktaki örnek yapılar.

**Yorum.** Eğim farkı büyük: VHH 33,4 Å²/kalıntı, VH 6,2 Å²/kalıntı (Δ = 27,2; q < 0,001) ve
her iki normalizasyon kontrolünden sonra da anlamlı kalıyor.

---

## 25 Ağustos, Salı — Figure 7: eski arayüzün antijenle etkileşimi

**Yapılan iş.** Aynı temas çekirdeği bu kez antikor–antijen sınırına uygulandı. Yapı başına
tek bir skaler üretiliyor: 25 pozisyonluk arayüzün kaçı antijenle temas ediyor (engagement
fraction). Mesafeler tek geçişte üç eşik için birden bayraklanıyor, böylece duyarlılık
analizi ek bir tarama gerektirmiyor:

```python
# src/qc/stage9_antigen_contacts.py
# gözlenmeyen pozisyonlar tabloda yok — sıfırlanmıyor (others are simply absent, not zeroed)
'contact_at_4A':   dist <= 4.0,
'contact_at_4p5A': dist <= 4.5,      # birincil eşik
'contact_at_5A':   dist <= 5.0,
engagement_fraction = n_engaged_positions / n_evaluable_positions
```

Görsel karşılaştırma için tek bir yapı seçmek yanıltıcı olurdu; bu yüzden **eşleştirilmiş
epitop** doğrulaması yazdım: aynı antijene bağlanan bir VHH ve bir VH–VL yapısı, antijen dizi
kimliği, temas kümesi Jaccard benzerliği ve antijen zincirlerinin süperpozisyonu (RMSD 2,34 Å,
N = 176 ortak CA atomu) ile eşleştirildi. İki panel aynı referans çerçevesinde çiziliyor,
yani görsel fark yalnızca antikor mimarisinden geliyor.

**Çıktı.** `Figure7.png` — (A) kohort dağılımı, (B/C) eşleştirilmiş IL4Rα kompleksleri ve
altında pozisyon-bazlı kalıntı şeritleri (kırmızı = antijenle temas eden arayüz pozisyonu).

**Yorum.** VHH'de medyan etkileşim oranı 0,36, VH'de 0,24; eşleştirilmiş örnekte de 10/25'e
karşı 7/25 — kohort istatistiği ile görsel örnek aynı yönde.

---

## 26 Ağustos, Çarşamba — Figure 9: yüzey fenotipi (PSH / PPC / PNC / elektrostatik)

**Yapılan iş.** Yük yamaları (PPC/PNC) dağılımlarının **sıfır-şişkin (zero-inflated)**
olduğunu fark ettim: yapıların büyük kısmında değer tam sıfır. Böyle bir dağılımda tek bir
ortalama testi hem "yama var mı" hem "ne kadar büyük" sorularını karıştırır. Bu yüzden metrik
iki ayrı teste bölündü ve asla tek kutu grafiğinde birleştirilmedi:

```python
# 1) yaygınlık: Fisher exact  → P(skor > 0) grup farkı
odds, p_prev = fisher_exact([[n_vhh_pos, n_vhh_neg], [n_vh_pos, n_vh_neg]])
# 2) koşullu büyüklük: yalnızca skor > 0 olan yapılar üzerinde Mann-Whitney U
u, p_mag = mannwhitneyu(vhh[vhh.PPC > 0].PPC, vh[vh.PPC > 0].PPC)
```

Elektrostatik katmanı için PDB2PQR + APBS harici araçları ayrı bir ortamda `subprocess` ile
çağrıldı, çıktı `.dx` potansiyel gridi okunup arayüz yüzeyinde örneklendi. Pozitif ve negatif
bileşenler **hiçbir zaman toplanmadı** (netlenmedi), çünkü net sıfır iki farklı fiziksel
durumu aynı sayıya indirger:

```python
# src/qc/stage8_apbs_pipeline.py
PDB2PQR_CMD = ["micromamba", "run", "-n", "vhh-elec", "pdb2pqr30", "--ff", "PARSE",
               "--apbs-input", ...]
APBS_CMD    = ["micromamba", "run", "-n", "vhh-elec", "apbs"]
if apbs_result.returncode != 0 or not dx_path.exists():
    return {'status': 'APBS_FAILURE', 'stderr': apbs_result.stderr[-2000:]}  # durum kodlu
```

**Çıktı.** `Figure9.png` — (A) PSH, (B) pozitif yama, (C) negatif yama, (D) APBS potansiyeli;
her satırda dağılım/oran grafiği + aynı örnek yapıların yüzey görselleri.

**Yorum.** VHH'de yük yaması bulundurma olasılığı belirgin şekilde yüksek (PPC %48,4 vs
%22,5; OR = 3,23), buna karşılık hidrofobik skor daha düşük — yani yüzey daha polar/yüklü bir
karaktere kaymış.

---

## 27 Ağustos, Perşembe — Figure 10: yerel dizi → özellik eşleme modülleri

**Yapılan iş.** O ana kadar altı farklı ölçüm (ıraksama, hidropati, yük, yüzey açıklığı, CDR3
gölgeleme, antijen etkileşimi) altı ayrı figürde duruyordu. Bu figürde hepsini **aynı IMGT
ekseninde**, dizi kompozisyonunun hemen altında sayısal değer olarak birleştirdim. Değerler
yeniden hesaplanmadı; dondurulmuş Stage-4 uzun tablosundan, aynı 533/1.429 temsilci kohortuna
süzülerek okundu:

```python
# src/qc/final_regional_physicochemistry.py
CHARGE_MAP = {'D': -1.0, 'E': -1.0, 'R': 1.0, 'K': 1.0, 'H': 0.1}   # source-verified, pH 7.4
FULL_EMPIRICAL_INTERFACE = {40, 42, 44, 49, 50, 52, 55, 66, 67, 68, 69, 103,
                            107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120}

df = long_df[long_df['INSTANCE'].isin(final_instances)].copy()
df = df[df['residue_present_in_sequence'] == True]   # koordinat-only kalıntılar dışarıda
df['charge'] = df['amino_acid'].map(CHARGE_MAP).fillna(0.0)
df['is_former_interface'] = df['imgt_position'].isin(FULL_EMPIRICAL_INTERFACE)
```

İki görselleştirme kuralı figürün yanıltıcı olmasını engelliyor. Birincisi, hangi
pozisyonların gösterileceği: her özellik için önceden ilan edilmiş "değere göre ilk 5" sırası
kullanıldı, liste sonradan elle değiştirilmedi. İkincisi, pencere birleştirme: her önemli
pozisyonun ±2 komşusu alınıyor, çakışan pencereler birleştiriliyor, ardışık olmayan bloklar
arasına görünür boşluk ve bölge etiketi konuyor — birbirinden uzak pozisyonlar asla yan
yanaymış gibi çizilmiyor. Metriğin tanımsız olduğu hücreler (örneğin CDR3 gölgelemesinin
CDR3'ün kendi içinde tanımsız olması) sıfırla doldurulmadı, taralı çizildi.

**Çıktı.** `Figure10.png` — altı modül; her modülde VHH logosu, VH logosu ve altında pozisyon
bazlı sayısal şerit.

**Yorum.** Aynı birkaç pozisyon (özellikle 42/49/50/52) birden fazla özellikte birlikte öne
çıkıyor; yani farklı ölçümler bağımsız değil, aynı bölgeye işaret ediyor.

---

## 28 Ağustos, Cuma — Figure 11: hotspot'ların 3B lokalizasyonu ve kapanış

**Yapılan iş.** Son figür, sayısal sonuçları tekrar yapı üzerine taşıyor. Buradaki asıl mesele
örnek yapı seçiminin **öznel olmaması**: her grup için çok değişkenli bir özellik uzayı
kuruldu, değişkenler z-skoruna çevrildi ve grubun kendi centroid'ine uzaklığı en küçük olan
yapılar seçildi. Beraberlik durumları için sabit bir sıralama tanımlı:

```python
# src/qc/final_freeze_3d_exemplars.py
def nearest(elig, metric_col, target, extra_tiebreak_cols=('resolution', 'SABDAB_ID')):
    e = elig.copy()
    e['dist'] = (e[metric_col] - target).abs()
    e = e.sort_values(['dist', *extra_tiebreak_cols])     # deterministik, rastgelelik yok
    return e.iloc[0]

# src/qc/stage12_raybould_alignment.py — ortak çerçeveye hizalama (SVD tabanlı)
def kabsch(P, Q):
    Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    rmsd = np.sqrt(np.mean(np.sum(((R @ Pc.T).T - Qc) ** 2, axis=1)))
    return R, Q.mean(axis=0) - R @ P.mean(axis=0), rmsd
```

Altı yapı ortak bir çerçeveye hizalandı ve 3B işaretleyicilerle alttaki dizi şeritleri
arasındaki renk eşlemesi otomatik olarak doğrulandı: **264 pozisyon-yapı kaydının tamamında
0 uyuşmazlık**.

**Çıktı.** `Figure11.png` — üstte 3 VHH, altta 3 VH yapısı; arayüz yarı saydam yüzey,
hotspot'lar altın küre, kanonik dörtlü magenta halka; her yapının altında kendi gerçek dizi
şeridi.

**Yorum ve staj kapanışı.** Üç haftada, veri alımından yayına hazır figürlere kadar uçtan uca
çalışan, her aşaması dosyaya dondurulmuş ve tekrar üretilebilir bir hat kuruldu. Öğrendiğim en
kalıcı şey, bir bilimsel pipeline'da kodun doğru çalışmasının yetmediği: eşiklerin tek bir
yerde tanımlanması, her seçimin deterministik bir kuralla yapılması, eksik verinin sıfırdan
ayrılması ve her ara sonucun denetim izi bırakması — sonucun güvenilirliğini belirleyen asıl
mühendislik kararları bunlar.
