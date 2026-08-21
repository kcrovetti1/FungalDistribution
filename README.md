# Global scale distribution of fungicide resistance

This project aims to better document the patterns of fungicide resistance at the global scale. We compared the unique accessions provided in Bédard et al. 2025. [^1], in which they built a global-scale dataset with over 64,00 entries. We manually curated their dataset filtering out entries with no accessions and without geographical information.

## Below is information regarding methods pertinent to the data in this repository. These methods are an excerpt from the app and manuscript submitted to...

**Title:**

# Global scale comparisons of fungicide resistance in clinical and agricultural settings.

Authors: Kincer Amburgey-Crovetti<sup>1\*</sup> and Teddy Garcia-Aroca<sup>1\*</sup>


<sup>1</sup> Department of Plant Pathology, University of Nebraska-Lincoln, Lincoln, NE 68503.


## Abstract
Fungicide resistance is a major issue at the global scale, with new groups of fungi across agricultural, non-agricultural, and clinical settings. We evaluated the global distribution of fungicide resistance across settings by manually curating the largest database for fungicide resistance called FungAMR. Additionally, we manually assembled a dataset for comparisons and to evaluate the accuracy of reporting across the globe. We encountered...
We also developed a [shiny app](https://kgacfel.shinyapps.io/FELFUNGCOMP/): https://kgacfel.shinyapps.io/FELFUNGCOMP/ with the goal of providing a resource that can be updated overtime.



| **CONTENTS**                                         |
| -----------------------------------------------------|
|												|
| 1. [INSTRUCTIONS](#instructions)														|
| 2. [DATA](#Data)                        |


# 1. INSTRUCTIONS 

In order to contribute, make changes, suggestions, or provide feedback to this repository, do the following:

1. [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) this repository from the top right.

<img width="479" alt="for_repo_image" src="https://docs.github.com/assets/cb-34352/mw-1440/images/help/repository/fork-button.webp">

2. Then, clone this repository to your laptop/desktop computer.

`
git clone https://github.com/kcrovetti1/FungalDistribution.git
`

3. Make changes, add, or edit current files and commit those changes to your copy of this 
repository:

| Command | Description |
| :--- | :------------------------------------- |
| `git status` | You should be able to see the changes you made (new files/folders) in red |
| `git add .` | To add all the changes to the current repository |
| `git commit -m "I changed/added x,y,z files"` | To commit those changes back to github |
| `git push` | To push the changes back to your forked repository |

Note: If you are having troubles when you try `git push`, follow [these](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) instructions to add your personal token in the command line.

4. Once you have pushed your changes, [submit a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests):

<img width="479" alt="pull-requests" src="https://docs.github.com/assets/cb-165497/mw-1440/images/help/pull_requests/merge-pull-request-options.webp">

Your changes will be reviewed by the repository owner (kcrovetti1) and accepted or denied depending on the accuracy/usefulness of the proposed changes.

# Data
This folder contains the metadata and manually curated datasets from FungAMR[^1] and the newly created dataset by the [Fungal Ecology Lab members](www.fungalecologylab.org): www.fungalecologylab.org.



[^1]:Bédard C, Pageau A, Fijarczyk A, Mendoza-Salido D, Alcañiz AJ, Després PC, Durand R, Plante S, Alexander EMM, Rouleau FD, Jordan DF, Jay A, Giguère M, Bernier M, Sharma J, Maroc L, Gervais NC, Menon ACT, Gagnon-Arsenault I, Bakker S, Rhodes J, Dufresne PJ, Bharat A, Sellam A, De Luca DG, Gerstein A, Shapiro RS, Quijada NM, Landry CR. FungAMR: a comprehensive database for investigating fungal mutations associated with antimicrobial resistance. Nat Microbiol. 2025 Sep;10(9):2338-2352. doi: 10.1038/s41564-025-02084-7. Epub 2025 Aug 11. PMID: 40790106.

