FLAGS = -Wall -Ofast
LINK = -llapack -lblas -llapacke -lm -ldl
SUBDIR = bin/
SRC = src/
OUT = out/
SUCCESS = @echo Build successful!

DIRS = $(shell find $(SRC) -type d)
SOURCES = $(foreach dir,$(DIRS),$(wildcard $(dir)/*.c))
OBJECTS = $(addprefix $(SUBDIR),$(SOURCES:.c=.o))

test: $(OBJECTS)
	@printf "\n\n$@: "
	@mkdir -p $(OUT)
	gcc $(OBJECTS) $(FLAGS) $(LINK) -o $(OUT)test
	$(SUCCESS)

clean:
	rm -r ./bin
	rm -r ./out

$(SUBDIR)%.o: %.c
	@printf "\n$(notdir $@): "
	@mkdir -p $(@D)
	gcc $(FLAGS) -c -o $@ $<
