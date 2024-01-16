import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Class for generating QC_plot from the FASTQ file
class QC_plot:

    def __init__(self, file_path):
        self.file_path = file_path

# Retrieve the position of the bases from the reads
    def get_base_positions(self):

        base_positions_list = []  # Empty-list to store the base_positions

        with open(self.file_path, 'r') as file:
            for i, line in enumerate(file):
                if i % 4 == 1:
                    base_position = list(range(1, len(line.strip()) + 1))
                    base_positions_list.extend(base_position)
        print({'Base_positions' : base_positions_list})
        return base_positions_list

# Retrieve the phred scores of the reads
    def get_phred_scores(self):

        phred_scores_list = []  # Empty-list to store the phred_scores

        with open(self.file_path, 'r') as file:
            for i, line in enumerate(file):
                if i % 4 == 3:
                    scores = line.strip()
                    for score in scores:
                        phred_33 = ord(score)-33  # Phred+33 encoding
                        phred_64 = ord(score)-64  # Phred+64 encoding

                        if ord(score)-33 >= 0 and ord(score)-33 <= 41:
                            phred_scores_list.append(phred_33)

                        elif ord(score)-64 >= 0 and ord(score)-64 <= 41:
                            phred_scores_list.append(phred_64)

                        else:
                            print("ERROR")

        return phred_scores_list


# File_path
result = QC_plot('/home/yamima/Python/FastQC/sample1.fastq')
bp = result.get_base_positions()
ps = result.get_phred_scores()




data = {'Base_position': bp, 'Phred_score': ps}
df = pd.DataFrame(data)
print(df)

# Box plot of per base sequence quality
sns.boxplot(data=df, x='Base_position', y='Phred_score',
            showfliers=False, color='yellow')

# Adding a line plot for the mean
sns.lineplot(data=df, x='Base_position', y='Phred_score',
             label='Mean', color='red', estimator=np.mean)

# Adding a line plot for the median
sns.lineplot(data=df, x='Base_position', y='Phred_score',
             label='Median', color='blue', estimator=np.median)

# Determining the step size for the x-axis ticks based on the maximun of the base_positions
if max(bp) < 40:
    step = 2
elif max(bp) < 100:
    step = 5
else:
    step = 10

# Setting up the x-axis ticks and y-axis ticks
plt.title('Per Base Sequence Quality')
plt.xlabel('position in bp')
plt.ylabel('phred_score')
plt.xticks(np.arange(0, max(bp)+1, step))
plt.yticks(np.arange(0, max(ps)+1, 2))
plt.savefig('Per_Base_Sequence_Quality.png')
plt.show()
