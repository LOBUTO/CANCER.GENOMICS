#r Functions

Function.scale.data.table<-function(d.t, col.protect=1:3, criteria=c(), scaling="SCALE"){
  #Scaling can be regular z-scoring ("SCALE") or 0-1 range ("ZERO")

  y.colnames <- colnames(d.t[,-col.protect,with=F])

  if (length(criteria)>0){
    all.criteria <- unique(d.t[[criteria]])

    all.y <- lapply(all.criteria, function(x)  {
      z <- d.t[d.t[[criteria]]==x,]
      z <- data.frame(z[,-col.protect,with=F])
      col.except <- apply(z, 2, sd)==0
      not.col.except <- apply(z, 2, sd)!=0

      if (scaling=="SCALE"){
        z <- cbind(z[,col.except], scale(z[,not.col.except]))
      } else if (scaling=="ZERO"){
        z <- cbind(z[,col.except], apply(z[,not.col.except], 2 ,function(zz) Function.range.0.1(zz)) )
      }

      return (z)
    })
    y <- do.call(rbind, all.y)

  } else {
    y<-data.frame(d.t[,-col.protect,with=F])

    if (scaling=="SCALE"){
      y<-scale(y)
    } else if (scaling=="ZERO"){
      y <- apply(y, 2, function(zz) Function.range.0.1(zz))
    }

  }

  x.colnames<-colnames(d.t[,col.protect,with=F])
  x<-cbind(d.t[,col.protect,with=F], y)
  setnames(x, c(x.colnames, y.colnames))

  return(x)
}
