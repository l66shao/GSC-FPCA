z=qnorm(.975)
for(v in c("controls","pooled")){d=read.csv(sprintf("latent_%s.csv",v)); cat(v,"\n"); for(k in split(d,sprintf("%04d %03d %02d",d$Nco,d$Ntr,d$m))){R=nrow(k);thd=k$ebar
 cvg=mean(k$plo<=0&k$phi>=0); frl=mean(k$flo<=thd&k$fhi>=thd); fth=mean(k$flo<=0&k$fhi>=0)
 th=mean(pnorm(z-k$B/k$pred_sd)-pnorm(-z-k$B/k$pred_sd))
 cat(sprintf(" %3d/%4d m=%2d reps=%3d | emp %.3f (se %.3f) | B %.3f | excess %.3f | sd %.3f predSD %.3f fitSD %.3f | pred/theta %.2f [formula %.2f] fit/real %.2f fit/theta %.2f | sigma2 %.2f\n",k$Ntr[1],k$Nco[1],k$m[1],R,mean(k$est),sd(k$est)/sqrt(R),mean(k$B),mean(k$est)-mean(k$B),sd(k$est),mean(k$pred_sd),mean(k$post_sd),cvg,th,frl,fth,mean(k$sigma2)))}}
