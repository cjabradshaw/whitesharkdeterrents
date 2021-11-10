# Response of white sharks exposed to newly developed personal shark deterrents

R code (by Corey Bradshaw @ Flinders University, Adelaide, Australia — <a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a>) to analyse deterrent trials for devices used to reduce the incidence of white shark (<em>Carcharodon carcharias</em>) attacks

Accompanies report:

C Huveneers<sup>1</sup>, S Whitmarsh<sup>1</sup>, M Thiele<sup>1</sup>, C May<sup>1</sup>, L Meyer<sup>1</sup>, A Fox<sup>2</sup>, CJA Bradshaw<sup>1</sup>

<sup>1</sup>College of Science and Engineering, Flinders University, Adelaide, South Australia
<sup>2</sup>Fox Shark Research Foundation, Adelaide, South Australia

Title: <a href="https://www.sharksmart.nsw.gov.au/__data/assets/pdf_file/0009/871785/Shark-response-to-personal-deterrents_Flinders.pdf"><em>Response of White Sharks Exposed to Newly Developed Personal Shark Deterrents</em></a>. Flinders University. Report to the NSW Department of Primary Industries.

*** NOTE *** This report has now been published in the peer-reviewed journal <em>PeerJ</em>: 

Huveneers, C, S Whitmarsh, M Thiele, L Meyer, A Fox, CJA Bradshaw. 2018. <a href="http://doi.org/10.7717/peerj.5554">Effectiveness of five personal shark-bite deterrents for surfers</a>. <em>PeerJ</em> 6: e5554 

## Abstract
The number of shark-human interactions and shark bites per capita has been increasing since the 1980s, leading to a rise in measures developed to mitigate the risk of shark bites. Yet many of the products commercially available for personal protection have not been scientifically tested, potentially providing an exaggerated sense of security to the people using them. We tested five personal shark deterrents developed for surfers (Shark Shield Pty Ltd [Ocean Guardian] Freedom+ Surf, Rpela, SharkBanz bracelet, SharkBanz surf leash, and Chillax Wax) by comparing the percentage of baits taken, distance to the bait, number of passes, and whether a shark reaction could be observed. We did a total of 297 successful trials at the Neptune Islands Group Marine Park in South Australia, during which 44 different white sharks (Carcharodon carcharias) interacted with the bait, making a total of 1413 passes. The effectiveness of the deterrents was variable, with the Freedom+ Surf affecting shark behaviour the most and reducing the percentage of bait taken from 96% (relative to the control board) to 40%. The mean distance of sharks to the board increased from 1.6 ± 0.1 m (control board) to 2.6 ± 0.1 m when the Freedom Surf+ was active. The other deterrents had limited or no measureable effect on white shark behavour. Based on our power analyses, the smallest effect size that could be reliably detected was ∼15%, which for the first time provides information about the effect size that a deterrent study like ours can reliably detect. Our study shows that deterrents based on similar principles—overwhelming a shark’s electroreceptors (the ampullae of Lorenzini) with electrical pulses—differ in their efficacy, reinforcing the need to test each product independently. Our results will allow private and government agencies and the public to make informed decisions about the use and suitability of these five products.


## Required R scripts
- <code>white shark deterrent trials Github.R</code> - base generalised linear mixed-effects models and accompanying power analyses
- <code>new_lmer_AIC_tables3.R</code> - functions to compare generalised linear mixed-effects models
- <code>r.squared.R</code> - functions to estimate R<sup>2</sup> for mixed-effects models

## Data files
NOTE: The following data files needed to run these analyses have now been uploaded. However, if you would like to use these files, contact <a href="mailto:charlie.huveneers@flinders.edu.au">Charlie Huveneers</a> or <a href="mailto:corey.bradshaw@flinders.edu.au">Corey Bradshaw</a>, explaining your intended uses for the data and code listed here, and to ask for permission. Thank you.<br>

- 'dist_data.csv'
- 'TF_data.csv'
- 'Binomial_data.csv'
- 'NoApp_data.csv'

###############################################################
