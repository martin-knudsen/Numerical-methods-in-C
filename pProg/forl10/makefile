blahblah : data log ; for f in $^; do echo -n $$f:; cat $$f; done
data log : makefile ; echo $@ > $@
clean    :          ; $(RM) data log