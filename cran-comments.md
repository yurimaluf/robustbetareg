## Resubmission
This is a resubmission.

## NOTES
> Please write TRUE and FALSE instead of T and F.
'T' and 'F' instead of TRUE and FALSE:
   man/plotenvelope.Rd:
     plotenvelope(
       object,
       type = c("sweighted2", "pearson", "weighted", "sweighted",
"sweighted.gamma",
         "sweighted2.gamma", "combined", "combined.projection"),
       conf = 0.95,
       n.sim = 100,
       PrgBar = T,
       control = robustbetareg.control(...),
       ...
     )

Done. 

> Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      methodsrobustbetareg.Rd: \value
      plot.robustbetareg.Rd: \value

Done.

> Some code lines in examples are commented out.
Please never do that. Ideally find toy examples that can be regularly
executed and checked. Lengthy examples (> 5 sec), can be wrapped in
\donttest{}.
Examples in comments in:
       plot.robustbetareg.Rd
       plotenvelope.Rd

Done. 

> Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.:  R/robustbetareg.R
If you're not familiar with the function, please check ?on.exit. This
function makes it possible to restore options before exiting a function
even if the function breaks. Therefore it needs to be called immediately
after the option change within a function.

Done. 

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies
There are currently no downstream dependencies for this package.
