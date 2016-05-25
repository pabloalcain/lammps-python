FOLDERS=$(wildcard */)
.PHONY: $(FOLDERS) all

all: $(FOLDERS)

$(FOLDERS):
	$(MAKE) -C $@

%:
	$(foreach DIR, $(FOLDERS), $(MAKE) -C $(DIR) $@;)
