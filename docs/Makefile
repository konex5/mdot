THEME_VERSION := v0.8.2
THEME := hugo-geekdoc
DESTINATION := target-doc
THEMEDIR := $(DESTINATION)/themes
CURRENTDIR := $(shell pwd)

.PHONY: doc
doc: doc-assets doc-build

.PHONY: doc-assets
doc-assets:
	mkdir -p $(THEMEDIR)/$(THEME)/ ; \
	curl -sSL "https://github.com/thegeeklab/$(THEME)/releases/download/${THEME_VERSION}/$(THEME).tar.gz" | tar -xz -C $(THEMEDIR)/$(THEME)/ --strip-components=1

.PHONY: doc-build
doc-build:
	cd $(DESTINATION); hugo --minify --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/ --baseURL docs/target-doc/public

.PHONY: clean
clean:
	rm -rf $(THEMEDIR) && \
	rm -rf $(DESTINATION)/public
serve:
	cd $(DESTINATION); hugo serve  --config=$(CURRENTDIR)/config.toml  --contentDir=$(CURRENTDIR)/content/
