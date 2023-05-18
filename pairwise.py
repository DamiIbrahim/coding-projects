#Citations - Shivani Rao & Arshia Kapil


import csv
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
from collections import deque


class CAM:
  """Find adhesion pairs given a scRNAseq adhesion molecule dataset."""

  def __init__(self, data, sig):
    self.data = data
    self.sig = sig
    self.pairs = []
    self.dfs = {0.1: 0, 0.05: 1, 0.025: 2, 0.01: 3, 0.005: 4, 0.0005: 5}
    self.critical_r = 0

  def correlation(self):
    """Generate correlation matrix."""

    df = pd.DataFrame(self.data)
    corr_matrix = df.corr(numeric_only=False)  # correlation matrix
    mask = np.zeros_like(corr_matrix)
    mask[np.triu_indices_from(mask)] = True
    sn.heatmap(corr_matrix, annot=True, mask=mask)
    return corr_matrix

  def critical(self):
    """Calculates and returns critical r."""

    N = len(list(self.data.values())[0])-2

    with open("rna_dataset - Critical.tsv", 'r') as file:
      next(file)
      csvreader = csv.reader(file)
      pearson_matrix = deque()

      for row in csvreader:
        item = row[0].split()
        item.pop(0)
        res = [eval(i) for i in item]
        pearson_matrix.append(res)
      pearson_matrix = np.array(pearson_matrix)
      self.critical_r = pearson_matrix[N-1][self.dfs.get(self.sig)]

    return self.critical_r

  def pairgen(self, corr_matrix, r):
    """Generates and returns adhesion pairs."""

    boolmatrix = corr_matrix.ge(self.critical_r)           #creates a boolean matrix depending on if matrix values less than or equal to critical r
    matches = boolmatrix[boolmatrix == True].stack().index.tolist()          #identifies the adhesion pairs
    matches = [match for match in matches if match[1] != match[0]]              #filters out CAMs that are equal to each other
    pairs = list({*map(tuple, map(sorted, matches))})          #filters out repeat CAMs
    self.pairs = [pair for pair in pairs if pair[0][0] == pair[1][0]]         #filter out PCDH:CDH pairs
    return self.pairs

  def plot(self):
    """Uses matplotlib to plot correlation matrix. """
    plt.show()

def main(inFile=None):
  """Calls the methods required for creating a correlation matrix and identifying valid cell adhesion pairs. """

  print("This tool generates cell adhesion pairs given an input of a scRNAseq tsv file in the following format:\n "
        "Gene ID, Gene name, Tissue, TPM, nTPM, pTPM ")
  print()

  try:
    sig = float(input("Enter significance value (alpha). Default 0.05. \nOptions: "
                      "0.1, 0.05, 0.025, 0.01, 0.005, 0.0005\n: "))
    if sig not in [0.1, 0.05, 0.025, 0.01, 0.005, 0.0005]:
      sig = 0.05
  except ValueError:
    sig = 0.05

  data = {}
  with open(inFile,'r') as file:
    next(file)
    csvreader = csv.reader(file)
    for row in csvreader:
      item = row[0].split()
      genetype = item[1]
      location = item[2]
      reads = float(item[-2])

      if genetype not in data.keys():
        data[genetype] = [reads]
      else:
        data[genetype].append(reads)

    thisCAM = CAM(data, sig)
    matrix = thisCAM.correlation()          #generate correlation matrix
    r = thisCAM.critical()              #determine critical r value
    ret_pairs = thisCAM.pairgen(matrix,r)         #generate and filter final adhesion pairs
    print("CAMs generated: ", ret_pairs)
    print(r)
    plot_matrix = thisCAM.plot()              #plot correlation matrix


if __name__ == "__main__":
    main("epithelial.tsv") #Name of tsv file for specific, user-inputted dataset