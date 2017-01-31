###################################
### Mathematica Code Parser     ###
### Author : Karol Nienałtowski ###
###################################

parser <- function(code){

rules <- rbind(
  c("r","reaction"),
  c(" ", "\\[RawSpace]"),
  c("!",  "\\[RawExclamation]"),
  c("\"",	"\\[RawDoubleQuote]"),
  c("#",	"\\[RawNumberSign]"),
  c("$", 	"\\[RawDollar]"),
  c("%", "\\[RawPercent]"),
  c("&", "\\[RawAmpersand]"),
  c("'", "\\[RawQuote]"),
  c("(", "\\[RawLeftParenthesis]"),
  c(")", "\\[RawRightParenthesis]"),
  c("*", "\\[RawStar]"),
  c("+", "\\[RawPlus]"),
  c(",", "\\[RawComma]"),
  c("-", "\\[RawDash]"),
  c(".", "\\[RawDot]"),
  c("/", "\\[RawSlash]"),
  c(":", "\\[RawColon]"),
  c(";", "\\[RawSemicolon]"),
  c("<", "\\[RawLess]"),
  c("=", "\\[RawEqual]"),
  c(">", "\\[RawGreater]"),
  c("?", "\\[RawQuestion]"),
  c("@", "\\[RawAt]"),
  c("[", "\\[RawLeftBracket]"),
  c("\\", "\\[RawBackslash]"),
  c("]", "\\[RawRightBracket]"),
  c("^", "\\[RawWedge]"),
  c("_", "\\[RawUnderscore]"),
  c("`", "\\[RawBackquote]"),
  c("{", "\\[RawLeftBrace]"),
  c("|", "\\[RawVerticalBar]"),
  c("}", "\\[RawRightBrace]"),
  c("~", "\\[RawTilde]"),
  c("_", "\\[UnderBrace]"),
  c("_", "\\[UnderBracket]"),
  c("_", "\\[UnderParenthesis]"),
  c("", "\" \\\""),
  c("", "\\\\"),
  c("", "\\\""),
  c(",", ", \""),
  c(",", ", "),
  c("", "\""),
  c(""," ")
)

parse <- function(code, i = 1){
  new_code <- gsub(rules[i,2], rules[i,1], code, fixed=TRUE)
  j <- i + 1
  ifelse(j < length(rules[,1]), 
         return(parse(new_code,j)),
         return(new_code))
}

result <- parse(code)
result <- unlist(strsplit(result, split=c(','), fixed = TRUE ))
result <- unlist(strsplit(result, split=c('{'), fixed = TRUE ))
result <- unlist(strsplit(result, split=c('}'), fixed = TRUE ))
result[result != ""]

}
