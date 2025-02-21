HUGO_OPTS := # Deployement: --baseURL https://nokx5.github.io/golden_cpp/
THEME_VERSION := v0.13.2
THEME_NAME := hugo-geekdoc
CURRENTDIR := $(shell pwd)
DESTINATION := $(CURRENTDIR)/target-doc
THEMEDIR := $(DESTINATION)/themes

.PHONY: doc
serve: clean doc-assets
	mkdir -p $(DESTINATION) && cd $(DESTINATION) && hugo serve  --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/ --themesDir=$(THEMEDIR) $(HUGO_OPTS) 

.PHONY: doc-assets
doc-assets:
	if ! [ -d $(THEMEDIR)/$(THEME_NAME) ]; then \
	  mkdir -p $(THEMEDIR)/$(THEME_NAME)/ ; \
	  curl -sSL "https://github.com/thegeeklab/$(THEME_NAME)/releases/download/${THEME_VERSION}/$(THEME_NAME).tar.gz" | tar -xz -C $(THEMEDIR)/$(THEME_NAME)/ --strip-components=1 ; \
	fi

.PHONY: doc-build
doc-build: doc-assets
	mkdir -p $(DESTINATION) && cd $(DESTINATION) && hugo --minify --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/ $(HUGO_OPTS) --themesDir=$(THEMEDIR)

.PHONY: clean
clean:
	rm -rf $(DESTINATION)/public

.PHONY: clean-all
clean-all:
	rm -rf $(DESTINATION)
