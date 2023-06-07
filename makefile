FLAGS = -Wall -O3
LINK = -lopenblas -lm -ldl -lscalapack-openmpi
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
	mpicc $(OBJECTS) $(FLAGS) $(LINK) -o $(OUT)test
	$(SUCCESS)

clean:
	rm -rf ./bin
	rm -rf ./out

$(SUBDIR)%.o: %.c
	@printf "\n$(notdir $@): "
	@mkdir -p $(@D)
	mpicc $(FLAGS) -c -o $@ $<
