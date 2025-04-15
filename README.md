# 🧬 Genome Mapinator

**Genome Mapinator** is a pipeline that automates genome annotation, feature visualization, and the creation of genome maps in PDF format. It combines tools like **Prodigal**, **AnnoMazing**, **Lovis4U**, and custom GFF manipulation scripts to produce a visual and annotated genome map from a folder of FASTA files.

---

## 📦 Features

- Runs gene prediction using **Prodigal**
- Annotates protein sequences with **Pfam** domains using **AnnoMazing**
- Edits and finalizes GFF files with associated sequence-region headers
- Generates annotated circular genome maps using **Lovis4U**
- Produces a clean final output with:
  - `Genome_map.pdf`
  - Annotated `.csv` files
  - Organized GFF and gene prediction files

---

## 🚀 Installation

Before using Genome Mapinator, make sure you have the following tools installed and accessible in your environment:

- Python 3.x with `pandas` and `numpy`
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [Lovis4U](https://github.com/moniruzzamanlab/lovis4u)
- [HMMER](http://hmmer.org/)

You will also need to download the Pfam database and place it in the hmms folder

---

## 🧰 Usage

```bash
python genome_mapinator.py -i /path/to/input_folder -o /path/to/output_folder
