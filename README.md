# FastSTest
YÜKSEK BOYUTLU GENOM VERİLERİNDE SIRALI ÖRÜNTÜLERİ SAHİP BAĞIMLI ÖZELLİKLERİN SEÇİMİ İÇİN S TESTİ UYARLAMASI
(Yüksek Lisans Tezi)
Deniz CEBELİ
GAZİ ÜNİVERSİTESİ
FEN BİLİMLERİ ENSTİTÜSÜ
Ağustos 2022
ÖZET
Yüksek boyutlu verilerde özellik seçimi makine öğrenmesindeki kritik adımlardan biridir.
Yüksek boyutlu veriler çok sayıda niteliğe karşın az sayıda gözlem içeren veri yapılarıdır.
Özellikle gen verilerine ilişkin çalışmalarda bu tarz verilerle çok sık karşılaşılmaktadır. Son
yıllarda makine öğrenmesi tekniklerinin yaygınlaşmasıyla genom çapında ilişkilendirme
çalışmaları (GWAS) artış göstermiştir. Bu tarz çalışmalarda tek nükleotid polimorfizm
(SNP) düzeyindeki artış ile marker değerlerindeki artış veya azalış örüntüleri tespit edilmeye
çalışılır. İstatistikte bu tarz örüntüler Jonckheere-Terpstra (JT), Terpstra-Magel (TM),
Ferdhiana-Terpstra-Magel (FTM), KTP, Modified JT ve S testi gibi sıralı alternatif
testleriyle incelenir. Ancak, yüksek boyutlu veriler için bu testlerin kullanımı hesaplama
zamanı bakımından ekonomik değildir. Bu nedenle, bu testlerin yüksek boyutlu veriler için
uyarlanması önem arz etmektedir. Bu çalışmada, aşırı çarpık dağılımlarda ve/veya
konveks/konkav alternatif hipotez durumlarında JT testine göre daha iyi sonuçlar veren S
istatistiğinin yüksek boyutlu veriler için uyarlanmış algoritması önerilmiştir. Elde edilen
sonuçlar S istatistiğinin yüksek boyutlu veriler için daha kullanışlı olduğunu göstermektedir.
fastS.R fonksiyonu: X ve Y değişkenleri için S istatistiğinin hesaplanması

fastSG.R fonksiyonu: Genom verisi için S istatistiğinin hesaplanması

generateData.R fonksiyonu: Simülasyon tasarımı
