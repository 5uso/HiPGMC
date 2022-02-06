FLAGS = -Wall -lblas -llapack -llapacke -Ofast
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
	gcc $(OBJECTS) $(FLAGS) -o $(OUT)test
	$(SUCCESS)

$(SUBDIR)%.o: %.c
	@printf "\n$(notdir $@): "
	@mkdir -p $(@D)
	gcc $(FLAGS) -c -o $@ $<
