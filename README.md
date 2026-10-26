# Reconstrução 3D da Bola — Método DLT

**Disciplina:** Biomecânica — EEFERP-USP  
**Tarefa 2**  
**Método:** Direct Linear Transformation (DLT) para reconstrução 3D  
**Referência:** Prof. Paulo R. P. Santiago — LaBioCoM-EEFERP-USP

---

## Descrição

Script `kvelball.py` que realiza a **reconstrução 3D da trajetória de uma bola chutada** a partir de imagens de **duas câmeras sincronizadas**, usando o método **DLT (Direct Linear Transformation)**.

O fluxo de trabalho é:

1. Digitalizar as coordenadas 2D da bola em cada câmera (ex.: com o software [Kinovea](https://www.kinovea.org/))
2. Digitalizar as coordenadas 2D dos pontos de calibração no calibrador volumétrico em cada câmera
3. Fornecer as coordenadas 3D reais dos pontos de calibração (`calibrador_ref.txt`)
4. Executar `kvelball.py` → obtém trajetória 3D, velocidade e ângulos de saída

---

## Estrutura dos arquivos

```
tarefa2/
├── kvelball.py            # Script principal de reconstrução 3D
├── calibrador_ref.txt     # Coordenadas 3D reais do calibrador (12 pontos)
├── requirements.txt       # Dependências Python
└── README.md              # Este arquivo

# Arquivos de entrada (gerados pelo Kinovea — não incluídos):
#   c1.txt       — coordenadas 2D da bola, câmera 1
#   c2.txt       — coordenadas 2D da bola, câmera 2
#   c1cal.txt    — coordenadas 2D do calibrador, câmera 1
#   c2cal.txt    — coordenadas 2D do calibrador, câmera 2
```

---

## Uso

### Instalação das dependências

```bash
pip install -r requirements.txt
```

### Execução

```bash
python kvelball.py c1.txt c2.txt c1cal.txt c2cal.txt calibrador_ref.txt <prefixo_saida>
```

**Exemplo:**
```bash
python kvelball.py c1.txt c2.txt c1cal.txt c2cal.txt calibrador_ref.txt resultado_chute
```

---

## Formato dos arquivos de entrada

### `c1.txt` / `c2.txt` (exportação Kinovea 0.9.4)

Arquivo com 4 linhas de cabeçalho, seguido de dados separados por espaço:

```
[cabeçalho linha 1]
[cabeçalho linha 2]
[cabeçalho linha 3]
[cabeçalho linha 4]
0  0  x0  y0
1  1  x1  y1
...
```

> As colunas utilizadas são as de índice 1 (frame) e 2 (x\_pixel, y\_pixel).

### `c1cal.txt` / `c2cal.txt` (coordenadas do calibrador no Kinovea)

Mesmo formato, sem linhas de cabeçalho (12 pontos):

```
x1  y1
x2  y2
...
x12 y12
```

### `calibrador_ref.txt` (coordenadas 3D reais)

12 pontos do calibrador volumétrico (em metros):

```
ponto  X      Y      Z
1      0      0      0
2      0      1.466  0
...
```

> O calibrador tem dimensões: 0,781 m × 1,466 m × 0,903 m (L × A × H)

---

## Arquivos de saída

| Arquivo                  | Conteúdo                                         |
|--------------------------|--------------------------------------------------|
| `<prefixo>.3d`           | Coordenadas 3D brutas (X, Y, Z em metros)        |
| `<prefixo>_filt.3d`      | Coordenadas 3D suavizadas (ajuste polinomial)    |
| `<prefixo>_result.txt`   | Velocidades instantâneas, azimute, elevação, média |

---

## Saídas numéricas (impressas no terminal)

- **Frame do impacto** detectado automaticamente
- **Velocidades instantâneas** (km/h) — entre frames após o impacto
- **Velocidade média** (km/h)
- **Azimute** (°) e **Elevação** (°) da trajetória de saída

---

## Método DLT — Resumo

O DLT mapeia coordenadas 3D do espaço (X, Y, Z) para coordenadas 2D na imagem (x, y) via 11 parâmetros L₁…L₁₁:

```
x = (L1·X + L2·Y + L3·Z + L4) / (L9·X + L10·Y + L11·Z + 1)
y = (L5·X + L6·Y + L7·Z + L8) / (L9·X + L10·Y + L11·Z + 1)
```

A calibração é feita por mínimos quadrados a partir de pontos de referência conhecidos. A reconstrução 3D usa os parâmetros de **duas câmeras** simultaneamente para calcular as coordenadas 3D.

---

## Parâmetros do experimento

| Parâmetro       | Valor          |
|----------------|----------------|
| Resolução Kinovea | 720 × 220 px |
| Frequência de aquisição | 120 Hz |
| Número de câmeras | 2 |
| Pontos de calibração | 12 |

---

## Dependências

- Python 3.x
- numpy
- pandas
- scipy
- matplotlib
