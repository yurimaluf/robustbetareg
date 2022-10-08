#'  Risk Managers Cost Effectiveness
#'
#'  The dataset are from a questionnaire that was sent to risk managers of large U.S.-based organizations.
#'
#' @format A data frame with 73 rows and 7 variables:
#' \describe{
#'   \item{FIRMCOST}{The measure of the firm’s risk management cost effectiveness, defined as total property and casualty premiums and uninsured losses as a percentage of total assets}
#'   \item{ASSUME}{Per occurrence retention amount as a percentage of total assets}
#'   \item{CAP}{Indicates that the firm owns a captive insurance company}
#'   \item{SIZELOG}{Logarithm of total assets}
#'   \item{INDCOST}{A measure of the firm’s industry risk}
#'   \item{CENTRAL}{A measure of the importance of the local managers in choosing the amount of risk to be retained}
#'   \item{SOPH}{A measure of the degree of importance in using analytical tools}
#'   ...
#' }
#'
#' @source \href{<https://instruction.bus.wisc.edu/jfrees/jfreesbooks/Regression\%20Modeling/BookWebDec2010/CSVData/RiskSurvey.csv>}{Schmit and Roth (1990)}
#'
#'@references \href{https://doi.org/10.2307/252842}{Schmit, Joan T., and Kendall Roth. “Cost Effectiveness of Risk Management Practices.” The Journal of Risk and Insurance, vol. 57, no. 3, 1990, pp. 455–70. JSTOR}
#'
#' @usage data("RiskManagerCost", package = "robustbetareg")
"RiskManagerCost"
