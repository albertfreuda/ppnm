LDLIBS = -lm

.PHONEY:default
default: out.cmdline.txt out.stdin.txt out.file.txt

out.cmdline.txt: read-cmdline input.txt
	./$< $$(cat input.txt) > out.cmdline.txt

out.stdin.txt: read-stdin input.txt
	./$< < input.txt > $@

out.file.txt: read-file input.txt
	./$< input.txt out.file.txt

input.txt: Makefile
	cat /dev/null > $@
	for x in $$(seq 0 0.5 5); do echo "$$x" >> $@; done

.PHONEY:clean
clean:
	$(RM) input.txt out* read-cmdline read-stdin read-file
