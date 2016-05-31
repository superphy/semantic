
/*
    Events
 */

(function() {
  var ContentHeader, CreateAccount, Factors, FactorsIntro, Foo, Foo_ooo, GeneResults, GeneResultsPanel, GeneSearchModel, GeneSearchPanel, GroupBrowse, Histogram, HistogramView, Home, LoginForm, Matrix, MatrixView, MetaData, Navbar, Page, Routes, Sidebar, SubmitSelectedView, SubmitView, Table, Table_2, Table_second, Test, Uploader, User, dragdrop, getEndpoint, getEndpoint2, submodule, url,
    indexOf = [].indexOf || function(item) { for (var i = 0, l = this.length; i < l; i++) { if (i in this && this[i] === item) return i; } return -1; },
    bind = function(fn, me){ return function(){ return fn.apply(me, arguments); }; },
    extend = function(child, parent) { for (var key in parent) { if (hasProp.call(parent, key)) child[key] = parent[key]; } function ctor() { this.constructor = child; } ctor.prototype = parent.prototype; child.prototype = new ctor(); child.__super__ = parent.prototype; return child; },
    hasProp = {}.hasOwnProperty;

  Routes = {};

  Routes.add = function(route, class_, args) {
    args = args || {};
    if (Routes.routes == null) {
      Routes.routes = {};
    }
    if (route.substring(0, 1) === '/') {
      route = route.substring(1);
    }
    return Routes.routes["/" + route] = {
      view: function() {
        return m.component(class_, args);
      }
    };
  };

  Page = (function() {
    function Page() {}

    Page.view = function(args) {
      args = args || {};
      return m('.', m.component(Navbar), args);
    };

    return Page;

  })();


  /*
  CLASS GeneSearchPanel
  
  Component class for the Gene Search Panel for one gene type
    in the VF and AMR gene feature.
  
  Args passed in:
      title: title of the panel
      type: string that describes type of gene (either "vf" or "amr")
      data: data model (global model that holds the selected genes & genomes)
      categories: object of categories for the genes
   */

  GeneSearchPanel = (function() {

    /*
    COMPONENT SelectedGenes
    
    Component that displays the selected genes in the search search panel.
    
    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        type: string that describes type of gene (either "vf" or "amr")
     */
    var CategorySelection, GeneTable, SearchSelect, SelectBar, SelectedGenes;

    function GeneSearchPanel() {}

    SelectedGenes = function(args) {
      var row;
      return m(".col-md-6 col-md-offset-3", [
        m(".selected-gene-list-wrapper", {
          id: args.type + "-selected-list"
        }, [
          m("fieldset", [
            m("span", ["Selected factors:"]), m("ul", {
              id: args.type + "-selected"
            }, [
              (function() {
                var k, len, ref, results1;
                ref = args.data.rows;
                results1 = [];
                for (k = 0, len = ref.length; k < len; k++) {
                  row = ref[k];
                  if (row.selected()) {
                    results1.push(m('li', {
                      "class": 'selected-gene-item'
                    }, [row["Gene_Name"]]));
                  }
                }
                return results1;
              })()
            ])
          ])
        ])
      ]);
    };


    /*
    COMPONENT GeneTable
    
    Component that displays the gene table.
    
    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        type: string that describes type of gene (either "vf" or "amr")
     */

    GeneTable = function(args) {
      return m(".col-md-6", [
        m(".gene-list-wrapper", [
          m("fieldset", [
            m("span.col-md-12", [
              "Select one or more " + args.type + " factors", m(".", m.component(Table, {
                data: args.data,
                checkbox: true
              }))
            ])
          ])
        ])
      ]);
    };


    /*
    COMPONENT SearchSelect
    
    Search bar for the gene table.
    
    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        type: string that describes type of gene (either "vf" or "amr")
     */

    SearchSelect = (function() {
      function SearchSelect() {}

      SearchSelect.controller = function(args) {
        var escapeRegExp;
        this.searchterm = m.prop("");
        escapeRegExp = function(str) {
          return str.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, "\\$&");
        };
        this.search = function(term) {
          var gene, k, len, ref, regex, row;
          self.searchterm = m.prop(term);
          regex = new RegExp(escapeRegExp(term), "i");
          console.log("Searching.....");
          ref = args.data.rows;
          for (k = 0, len = ref.length; k < len; k++) {
            row = ref[k];
            gene = row.Gene_Name;
            if (regex.test(gene)) {
              row.visible(true);
            } else {
              row.visible(false);
            }
          }
          return true;
        };
        this.select = function(all) {
          var k, len, ref, results1, row;
          ref = args.data.rows;
          results1 = [];
          for (k = 0, len = ref.length; k < len; k++) {
            row = ref[k];
            if (all === "true") {
              results1.push(row.selected(true));
            } else {
              results1.push(row.selected(false));
            }
          }
          return results1;
        };
        return this;
      };

      SearchSelect.view = function(ctrl, args) {
        return m(".gene-search-control-row", [
          m(".col-md-3", [
            m("input[type=text]", {
              id: args.type + "-autocomplete",
              "class": "form-control",
              placeholder: "Search gene in list",
              value: ctrl.searchterm(),
              onkeyup: m.withAttr("value", ctrl.search)
            })
          ]), m('.col-md-3', [
            m('.btn-group', [
              m('button', {
                "class": 'btn btn-link',
                checked: true,
                onclick: m.withAttr("checked", ctrl.select)
              }, "Select All"), m('button', {
                "class": 'btn btn-link',
                checked: false,
                onclick: m.withAttr("checked", ctrl.select)
              }, "Deselect All")
            ])
          ])
        ]);
      };

      return SearchSelect;

    })();


    /*
    COMPONENT CategorySelection
    
    Component for the category selection (whole section, not just one select bar)
    
    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        categories: object that has the categories and subcategories for the
                    particular set of genomes.
     */

    CategorySelection = (function() {
      function CategorySelection() {}

      CategorySelection.controller = function(args) {
        this.clear = function() {
          var k, len, ref, results1, row;
          ref = args.data.rows;
          results1 = [];
          for (k = 0, len = ref.length; k < len; k++) {
            row = ref[k];
            results1.push(row.visible(true));
          }
          return results1;
        };
        return this;
      };

      CategorySelection.view = function(ctrl, args) {
        var category;
        return m(".col-md-6", [
          m(".gene-category-wrapper", [
            m(".gene-category-intro", [
              m('span', "Select category to refine list of genes:"), m('button', {
                onclick: ctrl.clear
              }, "Clear")
            ]), (function() {
              var results1;
              results1 = [];
              for (category in args.categories) {
                results1.push([
                  m('.row', [m('.category-header col-xs-12', [category])]), m('.row', m('.col-xs-12', m.component(SelectBar, {
                    categories: args.categories,
                    category: category,
                    subcategories: args.categories[category],
                    data: args.data
                  })))
                ]);
              }
              return results1;
            })()
          ])
        ]);
      };

      return CategorySelection;

    })();


    /*
    COMPONENT SelectMult
    
    Component for each category search box. Based on the Select2 jQuery library
    
    Args:
        category: top category
        subcategories: list of subcategories
        data: gene search model
    
    WORK TO DO:
    - Fix problem with selecting multiple genes from different main categories
      Currently, if I pick subcategories from say "Adherence", then pick
      some others from "Autotransporter", the "Adherence" ones will disappear.
     */

    SelectBar = (function() {
      function SelectBar() {}

      SelectBar.controller = function(args) {
        this.allcategories = {};
        this.createCategories = function() {
          var cat, k, len, listCats, results1;
          listCats = Object.keys(args.categories);
          results1 = [];
          for (k = 0, len = listCats.length; k < len; k++) {
            cat = listCats[k];
            results1.push(self.allcategories[cat] = false);
          }
          return results1;
        };
        this.filterCategories = function(category, subcategories) {
          var cat, k, len, ref, results1, row, subcat;
          console.log("Filtering for ", category, " genes");
          ref = args.data.rows;
          results1 = [];
          for (k = 0, len = ref.length; k < len; k++) {
            row = ref[k];
            cat = row.Category;
            subcat = row.Sub_Category;
            if (!subcategories || (cat === category && indexOf.call(subcategories, subcat) >= 0)) {
              results1.push(row.visible(true));
            } else {
              results1.push(row.visible(false));
            }
          }
          return results1;
        };
        return this;
      };

      SelectBar.config = function(ctrl, args) {
        return function(element, isInitialized) {
          var attrs, el;
          ctrl.createCategories();
          el = $(element);
          if (!isInitialized) {
            attrs = {
              placeholder: "--Select a category--"
            };
            return el.select2(attrs).on("change", function(e) {
              var selected;
              m.startComputation();
              selected = [];
              if (el.val()) {
                selected = el.val();
              }
              ctrl.filterCategories(args.category, el.val());
              return m.endComputation();
            });
          }
        };
      };

      SelectBar.view = function(ctrl, args) {
        var attrs, index, subcat;
        return m("select", {
          "data-category-id": args.category,
          multiple: "multiple",
          "class": "form-control",
          config: SelectBar.config(ctrl, args)
        }, [
          (function() {
            var k, len, ref, results1;
            ref = args.subcategories;
            results1 = [];
            for (index = k = 0, len = ref.length; k < len; index = ++k) {
              subcat = ref[index];
              attrs = {
                id: "" + index,
                value: subcat,
                title: subcat
              };
              results1.push(m('option', attrs, subcat));
            }
            return results1;
          })()
        ]);
      };

      return SelectBar;

    })();

    GeneSearchPanel.controller = function(args) {
      return this;
    };

    GeneSearchPanel.view = function(args, ctrl) {
      return m(".panel panel-default", [
        m(".panel-heading", {
          id: ctrl.type + "-panel-header"
        }, [
          m("h4", {
            "class": "panel-title"
          }, m("a", {
            href: "#vf-form",
            config: m.route
          }, ctrl.title + " Form"))
        ]), m(".panel", {
          id: ctrl.type + "-panel"
        }, [
          m(".panel-body", [
            m(".row", [
              SelectedGenes({
                data: ctrl.data,
                type: ctrl.type
              })
            ]), m(".row", [
              m.component(SearchSelect, {
                type: ctrl.type,
                data: ctrl.data
              })
            ]), m(".row", [
              GeneTable({
                data: m.prop(ctrl.data),
                type: ctrl.type
              }), m.component(CategorySelection, {
                data: ctrl.data,
                categories: ctrl.categories
              })
            ])
          ])
        ])
      ]);
    };

    return GeneSearchPanel;

  })();

  HistogramView = (function() {
    function HistogramView() {
      this._create_histogram = bind(this._create_histogram, this);
      this.init = bind(this.init, this);
      return this;
    }

    HistogramView.prototype.init = function(searchResults, parentElem, elID) {
      var genes, genomes, view;
      genomes = Object.keys(searchResults);
      genes = Object.keys(searchResults[genomes[0]]);
      view = this;
      return function(el, isInitialized, ctx) {
        if (!isInitialized) {
          return view._create_histogram(searchResults, parentElem, elID);
        }
      };
    };

    HistogramView.prototype._create_histogram = function(searchResults, parentElem, elID) {
      var bins, margin;
      margin = {
        top: 40,
        right: 30,
        bottom: 40,
        left: 30
      };
      this.width = 300 - margin.left - margin.right;
      this.height = 250 - margin.top - margin.bottom;
      this.parentElem = parentElem;
      this.elID = elID;
      this.cssClass = "matrix_histogram";
      bins = [
        {
          'val': 0,
          'key': '0'
        }, {
          'val': 1,
          'key': '1'
        }, {
          'val': 2,
          'key': '2'
        }, {
          'val': 3,
          'key': '3'
        }, {
          'val': 4,
          'key': '4'
        }, {
          'val': 5,
          'key': '>=5'
        }
      ];
      this.x = d3.scale.ordinal().domain(bins.map(function(d) {
        return d.val;
      })).rangeRoundBands([0, this.width], .05);
      this.x2 = d3.scale.ordinal().domain(bins.map(function(d) {
        return d.key;
      })).rangeRoundBands([0, this.width], .05);
      this.xAxis = d3.svg.axis().scale(this.x2).orient("bottom");
      this.histogram = d3.layout.histogram().bins([0, 1, 2, 3, 4, 5, 6]);
      this.canvas = d3.select("#" + this.elID).append("svg").attr("width", this.width + margin.left + margin.right).attr("height", this.height + margin.top + margin.bottom).append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
      this.formatCount = d3.format(",.0f");
      this.canvas.append("g").attr("class", "x axis").attr("transform", "translate(0," + this.height + ")").call(this.xAxis).append("text").attr("dy", ".75em").attr("y", 23).attr("x", this.width / 2).attr("text-anchor", "middle").text('Number of Alleles');
      return this._getCounts(searchResults);
    };

    HistogramView.prototype._getCounts = function(searchResults) {
      var gene, genome, k, l, len, len1, ref, ref1, values;
      values = [];
      ref = Object.keys(searchResults);
      for (k = 0, len = ref.length; k < len; k++) {
        genome = ref[k];
        ref1 = Object.keys(searchResults[genome]);
        for (l = 0, len1 = ref1.length; l < len1; l++) {
          gene = ref1[l];
          values.push(searchResults[genome][gene]);
        }
      }
      console.log("values", values);
      return this.updateHistogram(values);
    };

    HistogramView.prototype.updateHistogram = function(values) {
      var histData, i, k, maxSteps, maxY, newBars, ref, steps, svgBars, yTop;
      histData = this.histogram(values);
      steps = [10, 50, 100, 200, 500, 800, 1000, 1200, 1500, 2000, 5000, 8000, 10000, 20000, 50000, 80000, 100000];
      maxSteps = steps.length;
      maxY = d3.max(histData, function(d) {
        return d.y;
      });
      yTop = NaN;
      for (i = k = 0, ref = maxSteps; k <= ref; i = k += 1) {
        if (maxY < steps[i]) {
          yTop = steps[i];
          break;
        }
      }
      this.y = d3.scale.linear().domain([0, yTop]).range([this.height, 0]);
      svgBars = this.canvas.selectAll("g.histobar").data(histData);
      svgBars.attr("transform", (function(_this) {
        return function(d) {
          return "translate(" + _this.x(d.x) + "," + _this.y(d.y) + ")";
        };
      })(this));
      svgBars.select("rect").attr("x", 0).attr("width", this.x.rangeBand()).attr("height", (function(_this) {
        return function(d) {
          return _this.height - _this.y(d.y);
        };
      })(this));
      svgBars.select("text").attr("dy", ".75em").attr("y", -14).attr("x", this.x.rangeBand() / 2).attr("text-anchor", "middle").text((function(_this) {
        return function(d) {
          if (d.y > 0) {
            return _this.formatCount(d.y);
          } else {
            return '';
          }
        };
      })(this));
      svgBars.exit().remove();
      newBars = svgBars.enter().append("g").attr("class", "histobar").attr("transform", (function(_this) {
        return function(d) {
          return "translate(" + _this.x(d.x) + "," + _this.y(d.y) + ")";
        };
      })(this));
      newBars.append("rect").attr("x", 0).attr("width", this.x.rangeBand()).attr("height", (function(_this) {
        return function(d) {
          return _this.height - _this.y(d.y);
        };
      })(this));
      newBars.append("text").attr("dy", ".75em").attr("y", -14).attr("x", this.x.rangeBand() / 2).attr("text-anchor", "middle").text((function(_this) {
        return function(d) {
          if (d.y > 0) {
            return _this.formatCount(d.y);
          } else {
            return '';
          }
        };
      })(this));
      return true;
    };

    return HistogramView;

  })();


  /*
  CLASS MatrixView
  
  Creates the matrix for the gene search results. 
  Uses D3 library.
  
  Args:
      None
   */

  MatrixView = (function() {
    function MatrixView() {
      this._compute_matrix = bind(this._compute_matrix, this);
      this._create_matrix = bind(this._create_matrix, this);
      this.init = bind(this.init, this);
      return this;
    }

    MatrixView.prototype.init = function(searchResults, parentElem, elID) {
      var genes, genomes, view;
      genomes = Object.keys(searchResults);
      genes = Object.keys(searchResults[genomes[0]]);
      view = this;
      return function(el, isInitialized, ctx) {
        if (!isInitialized) {
          return view._create_matrix(genomes, genes, searchResults, parentElem, elID);
        }
      };
    };

    MatrixView.prototype._create_matrix = function(genomes, genes, searchResults, parentElem, elID) {
      var dim, dropdown, dropdownDiv, num, self, wrap;
      console.log("Creating matrix...");
      self = this;
      this.cellWidth = 20;
      this.margin = {
        top: 150,
        right: 0,
        bottom: 0,
        left: 250
      };
      this.height = genomes.length * this.cellWidth;
      this.width = genes.length * this.cellWidth;
      dim = {
        w: this.width + this.margin.right + this.margin.left,
        h: this.height + this.margin.top + this.margin.bottom
      };
      this.parentElem = parentElem;
      this.elID = elID;
      this._compute_matrix(genomes, genes, searchResults);
      this.geneOrders = {
        name: d3.range(this.lengenes).sort((function(_this) {
          return function(a, b) {
            return d3.ascending(_this.geneNodes[a].name, _this.geneNodes[b].name);
          };
        })(this)),
        count: d3.range(this.lengenes).sort((function(_this) {
          return function(a, b) {
            return _this.geneNodes[b].count - _this.geneNodes[a].count;
          };
        })(this))
      };
      this.geneOrders['group'] = this.geneOrders['count'];
      this.orderType = 'name';
      this.z = d3.scale.linear().domain([0, 4]).clamp(true);
      this.x = d3.scale.ordinal().rangeBands([0, this.width]);
      this.cssClass = 'superphy-matrix';
      dropdownDiv = jQuery('<div class="matrixSort"><span>Order:</span> </div>').appendTo("#" + this.elID);
      dropdown = jQuery('<select name="matrix-sort">' + '<option value="name" selected="selected"> by Name</option>' + '<option value="count"> by Frequency</option>' + '<option value="group"> by Group</option>' + '</select>').appendTo(dropdownDiv);
      num = 0;
      dropdown.change(function() {
        var sortType;
        sortType = this.value;
        console.log("sortType", sortType);
        return self.viewAction(num, ['matrix-sort', sortType]);
      });
      wrap = d3.select("#" + this.elID).append("div").attr("class", "matrix-container").append("svg").attr("width", dim.w).attr("height", dim.h);
      this.canvas = wrap.append("g").attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");
      this.canvas.append("rect").attr("class", "matrixBackground").attr("width", this.width).attr("height", this.height).style("fill", "#f2f2f2");
      ({
        type: "matrix",
        elName: "genome_matrix",
        duration: 500
      });
      this.formatCount = d3.format(",.0f");
      return self.viewAction(num, ['matrix-sort', 'name']);
    };

    MatrixView.prototype._compute_matrix = function(genomes, genes, results) {
      var g, gObj, i, k, l, len, len1, len2, n, numAlleles, o, ref, results1;
      this.lengenomes = genomes.length;
      this.lengenes = genes.length;
      this.matrix = [];
      this.genomeNodes = [];
      i = 0;
      for (k = 0, len = genomes.length; k < len; k++) {
        g = genomes[k];
        gObj = {
          id: i,
          name: g,
          count: 0
        };
        this.genomeNodes.push(gObj);
        this.matrix[i] = d3.range(this.lengenes).map(function(j) {
          return {
            x: j,
            y: i,
            z: 0,
            i: null
          };
        });
        i++;
      }
      this.geneNodes = [];
      i = 0;
      for (l = 0, len1 = genes.length; l < len1; l++) {
        g = genes[l];
        gObj = {
          id: i,
          name: g,
          count: 0
        };
        this.geneNodes.push(gObj);
        i++;
      }
      i = 0;
      ref = this.genomeNodes;
      results1 = [];
      for (o = 0, len2 = ref.length; o < len2; o++) {
        g = ref[o];
        results1.push((function() {
          var len3, q, ref1, results2;
          ref1 = this.geneNodes;
          results2 = [];
          for (q = 0, len3 = ref1.length; q < len3; q++) {
            n = ref1[q];
            numAlleles = results[g.name][n.name];
            g.count += numAlleles;
            n.count += numAlleles;
            this.matrix[g.id][n.id].z = numAlleles;
            this.matrix[g.id][n.id].i = i;
            results2.push(i++);
          }
          return results2;
        }).call(this));
      }
      return results1;
    };

    MatrixView.prototype.viewAction = function(genomes, argArray) {
      var evt, ref;
      evt = argArray.shift();
      if (evt === 'matrix-sort') {
        this.orderType = argArray[0];
        if (ref = this.orderType, indexOf.call(Object.keys(this.geneOrders), ref) < 0) {
          throw new SuperphyError("Unrecognized order type: " + this.orderType + " in MatrixView viewAction method.");
        }
        return this._update(genomes);
      } else {
        throw new SuperphyError("Unrecognized event type: " + evt + " in MatrixView viewAction method.");
      }
    };

    MatrixView.prototype._update = function(genomes) {
      var genomesExit, newCols, newRows, svgGenes, svgGenomes, t1, that;
      t1 = new Date();
      this._sync(genomes);
      this.height = this.cellWidth * this.lengenomes;
      this.y = d3.scale.ordinal().rangeBands([0, this.height]);
      this.y.domain(this.genomeOrders[this.orderType]);
      this.x.domain(this.geneOrders[this.orderType]);
      this.canvas.selectAll(".matrixBackground").attr("height", this.height);
      svgGenomes = this.canvas.selectAll("g.matrixrow").data(this.genomeNodes, function(d) {
        return d.id;
      });
      svgGenomes.attr("class", (function(_this) {
        return function(d) {
          return _this._classList(d);
        };
      })(this)).select("text.matrixlabel").text(function(d) {
        return d.name;
      });
      svgGenomes.selectAll("g.matrixcell title").text(function(d) {
        return d.name + "title";
      });
      that = this;
      newRows = svgGenomes.enter().append("g").attr("class", (function(_this) {
        return function(d) {
          return _this._classList(d);
        };
      })(this)).attr("transform", function(d, i) {
        return "translate(0,0)";
      }).each(function(d) {
        return that._row(this, that.matrix[d.id], that.x, that.y, that.z);
      });
      newRows.append("line").attr("x2", this.width).style("fill", "white");
      newRows.append("text").attr("class", "matrixlabel").attr("x", -6).attr("y", this.y.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "end").text(function(d) {
        return d.name;
      }).on("click", function(d) {
        return jQuery('#dialog-matrix-row-select').data('row-id', d.genome).dialog('open');
      });
      svgGenes = this.canvas.selectAll("g.matrixcolumn").data(this.geneNodes, function(d) {
        return d.id;
      });
      svgGenes.selectAll("line").attr("x1", -this.height).style("color", "white");
      newCols = svgGenes.enter().append("g").attr("class", "matrixcolumn").attr("transform", function(d, i) {
        return "translate(" + 0 + ")rotate(-90)";
      });
      newCols.append("line").attr("x1", -this.height).style("color", "white");
      newCols.append("text").attr("class", "matrixlabel").attr("x", 6).attr("y", this.y.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "start").text(function(d) {
        return d.name;
      }).on("click", function(d) {
        return jQuery('#dialog-matrix-col-select').data('col-id', d.gene).dialog('open');
      });
      this._assumePositions();
      return genomesExit = svgGenomes.exit().transition().duration(this.duration).attr("transform", function(d) {
        return "translate(0,0)";
      }).remove();
    };

    MatrixView.prototype._sync = function(genomes) {
      var c, i, k, l, len, len1, n, ref, ref1;
      this.currNodes = [];
      this.currN = 0;
      ref = this.genomeNodes;
      for (k = 0, len = ref.length; k < len; k++) {
        n = ref[k];
        n.viewname = n.name;
        n.assignedGroup = 0;
        n.index = this.currN;
        this.currNodes.push(n);
        this.currN++;
        i = 0;
        ref1 = this.matrix[n.id];
        for (l = 0, len1 = ref1.length; l < len1; l++) {
          c = ref1[l];
          c.title = "genome title";
          i++;
        }
      }
      return this.genomeOrders = {
        name: d3.range(this.currN).sort((function(_this) {
          return function(a, b) {
            return d3.ascending(_this.currNodes[a].viewname, _this.currNodes[b].viewname);
          };
        })(this)),
        count: d3.range(this.currN).sort((function(_this) {
          return function(a, b) {
            return _this.currNodes[b].count - _this.currNodes[a].count;
          };
        })(this)),
        group: d3.range(this.currN).sort((function(_this) {
          return function(a, b) {
            var gdiff;
            gdiff = _this.currNodes[b].assignedGroup - _this.currNodes[a].assignedGroup;
            if (gdiff === 0) {
              return _this.currNodes[b].count - _this.currNodes[a].count;
            } else {
              return gdiff;
            }
          };
        })(this))
      };
    };

    MatrixView.prototype._row = function(svgRow, rowData, x, y, z) {
      var newCells, num, svgCells;
      svgCells = d3.select(svgRow).selectAll("matrixcell").data(rowData, function(d) {
        return d.i;
      });
      num = 0;
      newCells = svgCells.enter().append("g").attr("class", "matrixcell").attr("transform", (function(_this) {
        return function(d, i) {
          return "translate(" + _this.x(d.x) + ",0)";
        };
      })(this)).on("mouseover", (function(_this) {
        return function(p) {
          return _this._mouseover(p);
        };
      })(this)).on("mouseout", this._mouseout);
      newCells.append("rect").attr("x", 0).attr("width", x.rangeBand()).attr("height", y.rangeBand()).style("fill-opacity", function(d) {
        return z(d.z);
      });
      newCells.append("text").attr("dy", ".32em").attr("y", x.rangeBand() / 2).attr("x", x.rangeBand() / 2).attr("text-anchor", "middle").text((function(_this) {
        return function(d) {
          if (d.z > 0) {
            return _this.formatCount(d.z);
          } else {
            return '';
          }
        };
      })(this));
      return true;
    };

    MatrixView.prototype._assumePositions = function() {
      var that, transit;
      that = this;
      this.duration = 250;
      transit = this.canvas.transition().duration(this.duration);
      transit.selectAll(".matrixrow").attr("transform", function(d, i) {
        return "translate(0," + that.y(d.index) + ")";
      }).selectAll(".matrixcell").attr("x", function(d) {
        return that.x(d.x);
      }).attr("transform", function(d, i) {
        return "translate(" + that.x(d.x) + ",0)";
      });
      transit.selectAll(".matrixcolumn").attr("transform", function(d, i) {
        return "translate(" + that.x(i) + ")rotate(-90)";
      });
      return true;
    };

    MatrixView.prototype._mouseover = function(p) {
      d3.selectAll(".matrixrow text").classed("matrixActive", function(d, i) {
        return d.index === p.y;
      });
      return d3.selectAll(".matrixcolumn text").classed("matrixActive", function(d, i) {
        return i === p.x;
      });
    };

    MatrixView.prototype._mouseout = function() {
      return d3.selectAll("text").classed("matrixActive", false);
    };

    MatrixView.prototype._classList = function(d) {
      var clsList;
      clsList = ['matrixrow'];
      if (d.selected) {
        clsList.push("selectedRow");
      }
      if (d.assignedGroup != null) {
        clsList.push("groupedRow" + d.assignedGroup);
      }
      return clsList.join(' ');
    };

    return MatrixView;

  })();

  Navbar = (function() {
    function Navbar() {}

    Navbar.controller = function(args) {
      this.links = [
        {
          title: "Genome",
          url: "/Home.get()"
        }, {
          title: "Group Browse",
          url: "/gbrowse"
        }, {
          title: "Group Analyses",
          url: "/groups"
        }, {
          title: "VF and AMR",
          url: "/factors"
        }, {
          title: "Meta",
          url: "/meta"
        }
      ];
      return this;
    };

    Navbar.view = function(ctrl) {
      var link;
      return m(".container-fluid", m("nav.navbar navbar-default navbar-fixed-top", {
        role: 'navigation'
      }, m(".navbar-header", m("a.navbar-brand", {
        href: "/home",
        config: m.route
      }, "SuperPhy")), m("ul.nav navbar-nav", (function() {
        var k, len, ref, results1;
        ref = ctrl.links;
        results1 = [];
        for (k = 0, len = ref.length; k < len; k++) {
          link = ref[k];
          results1.push(m("li", m("a", {
            href: link.url,
            config: m.route
          }, link.title)));
        }
        return results1;
      })(), m("li.dropdown", m("a.dropdown-toggle", {
        href: "",
        role: "button",
        'data-toggle': "dropdown"
      }, "My Data", m("b.caret"))))));
    };

    return Navbar;

  })();


  /*
  CLASS Sidebar
  
  Component class for sidebar in VF and AMR feature
  
  Args passed in:
      None yet
  
  To do: Add components into Sidebar
   */

  Sidebar = (function() {
    function Sidebar() {}

    Sidebar.controller = function(args) {
      return this;
    };

    Sidebar.view = function(ctrl) {
      return m('.', {
        id: 'sidebar-wrapper'
      }, [
        m('.sidebar', [
          m('#sidebar-expand-collapse', [
            m('button.btn btn-default[type=button]', {
              id: 'sidebar-expand'
            }, ">>"), m('button.btn btn-default[type=button]', {
              id: 'sidebar-collapse'
            }, "<<")
          ])
        ])
      ]);
    };

    return Sidebar;

  })();


  /*
  CLASS Table
  
  Component class for table views
  
  Args passed in:
      data: Data to populate the table
      checkbox: boolean for checkbox option (for VF & AMR)
  
  To fix: make the VF & AMR checkbox portion look nicer?
   */

  Table = (function() {
    var state;

    function Table() {}

    state = {
      pageY: 0,
      pageHeight: window.innerHeight
    };

    window.addEventListener("scroll", function(e) {
      state.pageY = Math.max(e.pageY || window.pageYOffset, 0);
      state.pageHeight = window.innerHeight;
      return m.redraw();
    });

    Table.controller = function(args) {
      this.sort_table = function(list, attribute) {
        if (attribute == null) {
          attribute = 'data-sort-by';
        }
        return {
          onclick: function(e) {
            var first, item;
            item = e.target.getAttribute(attribute);
            if (item) {
              first = list[0];
              list.sort(function(a, b) {
                if (isNaN(parseFloat(a[item] * 1))) {
                  if (isNaN(parseFloat(b[item] * 1))) {
                    if (a[item] > b[item]) {
                      return 1;
                    } else if (b[item] > a[item]) {
                      return -1;
                    } else {
                      return 0;
                    }
                  } else {
                    return -1;
                  }
                } else if (isNaN(parseFloat(b[item] * 1))) {
                  return 1;
                } else if (a[item] * 1 < b[item] * 1) {
                  return 1;
                } else if (b[item] * 1 < a[item] * 1) {
                  return -1;
                } else {
                  return 0;
                }
              });
              if (first === list[0]) {
                list.reverse();
              }
            }
          }
        };
      };
      return this;
    };

    Table.view = function(ctrl, args) {
      var begin, end, header, headers, list, offset, pageY, row, rows, x;
      rows = args.data().rows;
      headers = args.data().headers;
      pageY = state.pageY;
      begin = pageY / 60 | 0;
      end = begin + (state.pageHeight / 60 | 0 + 10);
      offset = pageY % 60;
      return m(".Occlusion", {
        style: {
          height: args.data().rows.length * 46 + "px",
          position: "relative",
          top: -offset + "px"
        }
      }, [
        m("table", {
          style: {
            top: state.pageY + "px"
          }
        }, [
          m("tr", [
            (function() {
              var k, len, ref, results1;
              ref = args.data().headers;
              results1 = [];
              for (k = 0, len = ref.length; k < len; k++) {
                header = ref[k];
                results1.push(m('th[data-sort-by=' + header + ']', ctrl.sort_table(list = args.data().rows), [header]));
              }
              return results1;
            })()
          ]), (function() {
            var k, len, ref, results1;
            ref = args.data().rows.slice(begin, +end + 1 || 9e9);
            results1 = [];
            for (x = k = 0, len = ref.length; k < len; x = ++k) {
              row = ref[x];
              results1.push(m('tr', [
                (function() {
                  var l, len1, len2, o, ref1, results2, results3;
                  if (args.checkbox) {
                    results2 = [];
                    for (l = 0, len1 = headers.length; l < len1; l++) {
                      header = headers[l];
                      if (header === "Gene") {
                        results2.push(m('td', {
                          "class": 'gene_table_item'
                        }, [
                          m('.checkbox', [
                            m('label', [
                              m('input[type=checkbox]', {
                                "class": 'checkbox gene-table-checkbox gene-search-select',
                                checked: row.selected(),
                                onclick: m.withAttr("checked", row.selected)
                              }), row[header]
                            ])
                          ])
                        ]));
                      } else {
                        results2.push(m('td', {
                          "class": 'gene_table_item'
                        }, [m('label', [row[header]])]));
                      }
                    }
                    return results2;
                  } else {
                    ref1 = args.data().headers;
                    results3 = [];
                    for (o = 0, len2 = ref1.length; o < len2; o++) {
                      header = ref1[o];
                      results3.push(m('td', [row[header]]));
                    }
                    return results3;
                  }
                })()
              ]));
            }
            return results1;
          })()
        ])
      ]);
    };

    return Table;

  })();

  CreateAccount = (function(superClass) {
    var FormGroup;

    extend(CreateAccount, superClass);

    function CreateAccount() {
      return CreateAccount.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/SignUp', CreateAccount);

    FormGroup = (function() {
      function FormGroup() {}

      FormGroup.controller = function(args) {
        args = args || {};
        this.placeholder = args.placeholder || "";
        this.type = args.type || "text";
        this.value = args.value || m.prop("");
        this.help_block = args.help_block || null;
        return this;
      };

      FormGroup.view = function(ctrl) {
        return m('.form-group', m('label.col-sm-2 control-label', ctrl.placeholder), m('.col-sm-8', m('input.form-control', {
          type: ctrl.type,
          oninput: m.withAttr('value', ctrl.value),
          value: ctrl.value(),
          placeholder: ctrl.placeholder
        }), ctrl.help_block ? m('span.help-block', "(" + ctrl.help_block + ")") : void 0));
      };

      return FormGroup;

    })();

    CreateAccount.controller = function(args) {
      this.user = new User();
      this.submit = function() {
        return User.sign_up(this.user);
      };
      return this;
    };

    CreateAccount.view = function(ctrl) {
      return CreateAccount.__super__.constructor.view.call(this, m('.container', m('.content-header', m('h1', m('span.title_part1', "CREATE"), m('span.title_part2', "ACCOUNT"))), m('.panel panel-default', m('panel-body', m('form.form-horizontal', m('.form-group', m('.col-sm-offset-2 col-sm-8', m('p', "Already registered or want to make changes to your existing account?", m('a', {
        tabindex: "9",
        href: "edit_account"
      }, "Sign in")), m('p', "To create a new account, fill in all fields below."))), m.component(FormGroup, {
        value: ctrl.user.username,
        type: "text",
        placeholder: "Username",
        help_block: "alphanumeric"
      }), m.component(FormGroup, {
        value: ctrl.user.password,
        type: "Password",
        placeholder: "Password",
        help_block: "6-10 characters"
      }), m.component(FormGroup, {
        value: ctrl.user.password2,
        type: "Password",
        placeholder: "Re-enter Password"
      }), m.component(FormGroup, {
        value: ctrl.user.first_name,
        type: "text",
        placeholder: "First Name"
      }), m.component(FormGroup, {
        value: ctrl.user.last_name,
        id: "inputLastname",
        type: "text",
        placeholder: "Last Name"
      }), m.component(FormGroup, {
        value: ctrl.user.email,
        id: "inputEmail",
        type: "text",
        placeholder: "Email"
      }), m('.form-group', m('.col-sm-offset-2 col-sm-8', m("button.btn btn-primary[type=button]", {
        onclick: function() {
          return ctrl.submit();
        }
      }, "Create Account"))))))));
    };

    return CreateAccount;

  })(Page);

  LoginForm = (function(superClass) {
    extend(LoginForm, superClass);

    function LoginForm() {
      return LoginForm.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/Login', LoginForm);

    LoginForm.controller = function(args) {
      this.user = new User();
      this.submit = function() {
        return User.log_in(this.user);
      };
      return this;
    };

    LoginForm.view = function(ctrl) {
      return LoginForm.__super__.constructor.view.call(this, m('.container', m('.panel panel-default', m('panel-body', m('form.form-horizontal', m('h3', "Log in"), m('input[type=text]', {
        oninput: m.withAttr('value', ctrl.user.username),
        value: ctrl.user.username(),
        placeholder: "Username"
      }), m('input[type=Password]', {
        oninput: m.withAttr('value', ctrl.user.password),
        value: ctrl.user.password(),
        placeholder: "Password"
      }), m("button.btn btn-primary[type=button]", {
        onclick: function() {
          return ctrl.submit();
        }
      }, "Login"))))));
    };

    return LoginForm;

  })(Page);


  /*
  CLASS Factors
  
  Page component for the Factors (VF and AMR) page
   */

  Factors = (function(superClass) {
    extend(Factors, superClass);

    function Factors() {
      return Factors.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/factors', Factors);

    Factors.controller = function(args) {
      this.active = m.prop("genes");
      this.tabCtrl = new mc.Tabs.controller('genes');
      this.model = GeneSearchModel;
      this.vfcategories = this.model.getCategories('vf');
      this.amrcategories = this.model.getCategories('amr');
      return this;
    };

    Factors.view = function(ctrl) {
      var renderGeneForm, renderGenomeSearch, renderSubmit, renderTabContents;
      renderTabContents = function(ctrl) {
        var activeTab;
        activeTab = ctrl.tabCtrl.activeTabName;
        switch (activeTab()) {
          case 'genes':
            return renderGeneForm(ctrl);
          case 'genomes':
            return renderGenomeSearch(ctrl);
          case 'submit':
            return renderSubmit(ctrl);
          default:
            return m('p', "Not a valid tab name.");
        }
      };
      renderGeneForm = function(ctrl) {
        return m('.tab-content', [
          m('.tab-pane active', {
            id: 'gene-search-querygenes'
          }, [
            m('.panel-group genes-search', {
              id: 'accordian'
            }, [
              m.component(GeneSearchPanel, {
                title: "Virulence Factor",
                type: "vf",
                data: ctrl.model.vfList,
                categories: ctrl.vfcategories()
              }), m.component(GeneSearchPanel, {
                title: "Antimicrobial Resistance",
                type: "amr",
                data: ctrl.model.amrList,
                categories: ctrl.vfcategories()
              })
            ])
          ])
        ]);
      };
      renderGenomeSearch = function(ctrl) {
        return m('.tab-content', [
          m('.tab-pane active', {
            id: 'gene-search-querygenes'
          }, [
            m('.panel-group genes-search', {
              id: 'accordian'
            }, ["Testing Genome Search (you can change the genomes manually in the GeneSearchModel)"])
          ])
        ]);
      };
      renderSubmit = function(ctrl) {
        return m('.tab-content', [
          m('.tab-pane active', {
            id: 'gene-search-querygenes'
          }, [
            m('.panel-group genes-search', {
              id: 'accordian'
            }, [
              m.component(SubmitView, {
                data: ctrl.model
              })
            ])
          ])
        ]);
      };
      return Factors.__super__.constructor.view.call(this, m('.', {
        id: "wrapper"
      }, [
        m('.', {
          id: 'page-content-wrapper'
        }, [
          m('.', {
            id: 'page-content -inset'
          }, [
            m('.container-fluid', [
              m.component(FactorsIntro), m('.container', [
                mc.Tabs.view(ctrl.tabCtrl, {
                  flavor: 'bs/nav-tabs',
                  tabs: [
                    {
                      name: 'genes',
                      label: 'Select Genes'
                    }, {
                      name: 'genomes',
                      label: 'Select Genomes'
                    }, {
                      name: 'submit',
                      label: 'Submit Query'
                    }
                  ]
                }), renderTabContents(ctrl)
              ])
            ])
          ])
        ])
      ]));
    };

    return Factors;

  })(Page);


  /*
  COMPONENT FactorsIntro
  
  Introduction for gene indentification feature
  
  Args:
      none
   */

  FactorsIntro = {
    view: function(ctrl, args) {
      return m('.intro', [
        m('.row', [
          m('.col-xs-8', [m('.content-header', [m('h1', [m('span.title_part1', 'VIRULENCE & AMR '), m('span.title_part2', 'GENES')])])]), m('.col-xs-4', [
            m('button.btn btn-danger btn-lg', {
              id: "intro-button",
              type: "button"
            }, ["INTRODUCTION"])
          ])
        ]), m('p', "Search for the presence or absence of virulence factor genes or antimicrobial resistance genes in genomes of interest. Detailed information on individual virulence factor or antimicrobial resistance genes can be retrieved by clicking on the individual genes.")
      ]);
    }
  };


  /*
  COMPONENT SubmitView
  
  View for the third tab on the page 
  (This could be moved to another file?)
  
  Args:
      data: gene search model
   */

  SubmitView = {
    controller: function(args) {
      this.submit = function() {
        return args.data.submit();
      };
      this.reset = function() {
        return args.data.reset();
      };
      return this;
    },
    view: function(ctrl, args) {
      args.data.getSelectedVF();
      args.data.getSelectedAMR();
      return m('.tab-content', [
        m('.tab-pane active', {
          id: 'gene-search-submit'
        }, [
          m.component(SubmitSelectedView, {
            selected: args.data.selectedVF,
            title: "Virulence Factor"
          }), m.component(SubmitSelectedView, {
            selected: args.data.selectedAMR,
            title: "Antimicrobial Resistance"
          }), m.component(SubmitSelectedView, {
            selected: args.data.selectedGenomes,
            title: "Genome"
          }), m('.row', [
            m('.gene-search-next-wrapper', {
              id: 'query-gene-form'
            }, [
              m('button', {
                "class": 'btn btn-success',
                type: 'submit',
                value: 'Submit',
                onclick: ctrl.submit
              }, "Submit"), m('button', {
                "class": 'btn btn-danger',
                type: 'reset',
                value: 'Reset',
                onclick: ctrl.reset
              }, "Reset")
            ])
          ])
        ])
      ]);
    }
  };


  /*
  COMPONENT SubmitSelectedView
  
  Component that confirms number of gene selected
  
  Args:
      selected: list of selected genes or genomes
      title: type of selected genes or genomes
   */

  SubmitSelectedView = {
    view: function(ctrl, args) {
      return m('.row', [
        m('.col-md-4 col-md-offset-1', [
          m('.panel panel-default', [
            m('.panel-body', {
              id: 'vf-selected-count'
            }, [args.selected.length, " ", args.title.toLowerCase(), args.selected.length !== 1 ? "s" : void 0, " selected"])
          ])
        ])
      ]);
    }
  };


  /*
  CLASS GeneResults
  
  Page Component for the gene results page.
   */

  GeneResults = (function(superClass) {
    extend(GeneResults, superClass);

    function GeneResults() {
      return GeneResults.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/results', GeneResults);

    GeneResults.controller = function(args) {
      this.results = GeneSearchModel;
      return this;
    };

    GeneResults.view = function(ctrl, args) {
      var data, testdata;
      testdata = {
        genomes: ["JHNV00000000", "ANVW00000000"],
        vfs: ["saa", "papC", "ompA", "hylA", "QnrS7"],
        amrs: ["QnrS7"]
      };
      data = {
        vfresults: GeneSearchModel.vfresults,
        amrresults: GeneSearchModel.amrresults
      };
      return GeneResults.__super__.constructor.view.call(this, m('.', {
        id: 'wrapper'
      }, m.component(Sidebar), m('.', {
        id: "page-content-wrapper"
      }, m('.page-content inset', m('.container-fluid', m.component(ContentHeader, {
        title: "Virulence Factor and AMR Results"
      }), m('.', {
        id: 'results'
      }, [
        m.component(GeneResultsPanel, {
          id: "vf_results",
          type: "Virulence Factor",
          numSelected: GeneSearchModel.selectedVF.length,
          matrix: m.component(Matrix, {
            matrixview: new MatrixView(),
            results: data.vfresults,
            parentEl: "vf_result_matrix",
            elID: "genome_matrix1"
          }),
          histogram: m.component(Histogram, {
            histview: new HistogramView(),
            results: data.vfresults,
            parentEl: "vf_result_histogram",
            elID: "matrix_ticker1"
          })
        }), m.component(GeneResultsPanel, {
          id: "amr_results",
          type: "Antimicrobrial Resistance",
          numSelected: GeneSearchModel.selectedAMR.length,
          matrix: m.component(Matrix, {
            matrixview: new MatrixView(),
            results: data.amrresults,
            parentEl: "amr_result_matrix",
            elID: "genome_matrix2"
          }),
          histogram: m.component(Histogram, {
            histview: new HistogramView(),
            results: data.amrresults,
            parentEl: "amr_result_histogram",
            elID: "matrix_ticker2"
          })
        })
      ]))))));
    };

    return GeneResults;

  })(Page);


  /*
  COMPONENT GeneResultsPanel
  
  Wrapper component for results of a set of genes.
  
  Args:
      id: ID name for div element
      type: type of gene
      numSelected: number of selected genes
      matrix: MatrixView component
      histogram: Histogram component
   */

  GeneResultsPanel = {
    view: function(ctrl, args) {
      if (args.numSelected > 0) {
        return m('div', {
          id: args.id
        }, [m('hr'), m('h4', "Detected " + args.type + " Alleles"), args.matrix, args.histogram]);
      } else {
        return m('div', '');
      }
    }
  };


  /*
  COMPONENT ContentHeader
  
  Intro information for gene results
  
  Args:
      title: Title of the page
   */

  ContentHeader = {
    view: function(ctrl, args) {
      return m('.row', m('.col-xs-8', m('.content-header', m('h1', args.title))), m('.col-xs-4', m('button', {
        id: "intro-button",
        "class": "btn btn-danger-lg",
        type: "button"
      }, "INTRODUCTION")));
    }
  };


  /*
  COMPONENT Matrix
  
  Initializer component for the Matrix
  
  Args:
      matrixview: MatrixView component
      results: results object
      parentEl: parent element's ID
      elID: matrix's element ID
   */

  Matrix = {
    controller: function(args) {
      this.results = args.results();
      return this;
    },
    view: function(ctrl, args) {
      return m('div', {
        id: args.parentEl
      }, m('.superphy-matrix', {
        id: args.elID,
        config: args.matrixview.init(ctrl.results, args.parentEl, args.elID)
      }));
    }
  };


  /*
  COMPONENT Histogram
  
  Initializer component for the Histogram
  
  Args:
      histview: Histrogram component
      results: results object
      parentEl: parent element's ID
      elID: histogram's element ID
   */

  Histogram = {
    controller: function(args) {
      this.results = args.results();
      return this;
    },
    view: function(ctrl, args) {
      return m('.row histogram-row', m('.col-md-4 histogram-description'), m('span', "Blurb goes here"), m('.col-md-8 histogram-container', {
        id: args.parentEl
      }, m('.superphy-histogram', {
        id: args.elID,
        config: args.histview.init(ctrl.results, args.parentEl, args.elID)
      })));
    }
  };

  GroupBrowse = (function(superClass) {
    extend(GroupBrowse, superClass);

    function GroupBrowse() {
      return GroupBrowse.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/gbrowse', GroupBrowse);

    GroupBrowse.controller = function(args) {
      var url;
      this.data = getEndpoint2(url = "data/meta");
      return this;
    };

    GroupBrowse.view = function(ctrl) {
      return GroupBrowse.__super__.constructor.view.call(this, m(".", m.component(Table, {
        data: ctrl.data
      })));
    };

    return GroupBrowse;

  })(Page);


  /*
  CLASS Home
  
  Page component for the main page
   */

  Home = (function(superClass) {
    extend(Home, superClass);

    function Home() {
      return Home.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/', Home);

    Routes.add('/home', Home);

    Home.controller = function(args) {
      args = args || {};
      this.debug = args.debugging || true;
      return this;
    };

    Home.view = function(ctrl) {
      return Home.__super__.constructor.view.call(this, m(".container", {
        id: 'Home.get()-beta'
      }, ctrl.debug === true ? m("input[type=button]", {
        onclick: function() {
          localStorage.clear();
          return console.log("Cache Cleared!");
        },
        value: "Clear Cache"
      }) : void 0, m(".row", m(".well center-block", m("p.text-center", m("span.text-info beta-release", "Beta Release"), " Some features are still under development and may not be fully functional."))), m(".row superphy-image", m("img", {
        src: 'images/superphy_logo_with_title.png'
      }), m("p.superphy-image", 'NEXT-LEVEL PHYLOGENETIC AND EPIDEMIOLOGICAL ANALYSIS OF PATHOGENS')), m(".row well", m("p", 'A user-friendly, integrated platform for the predictive genomic analyses of ', m("em", 'Escherichia coli')), m("p", 'Provides near real-time analyses of thousands of genomes using novel computational approaches.'), m("p", 'Generates results that are understandable and useful to a wide community of users.')), m(".row", m("button.btn btn-danger btn-lg center-block", "INTRODUCTION")), m(".row text-center", m("p", 'For more information or to discuss possible collaboration, please contact:', m("p", m("ul", m("li", 'Dr. Vic Gannon: vic.gannon@phac-aspc.gc.ca')))))));
    };

    return Home;

  })(Page);


  /*
  Table_Second.coffee
      A trial table to test a different pagination style.
  
  Args passed in:
      data: Data to populate the table
   */

  Table_2 = (function() {
    var state;

    function Table_2() {}

    state = {
      pageY: 0,
      pageHeight: window.innerHeight
    };

    window.addEventListener("scroll", function(e) {
      state.pageY = Math.max(e.pageY || window.pageYOffset, 0);
      state.pageHeight = window.innerHeight;
      return m.redraw();
    });

    Table_2.controller = function(args) {
      this.sort_table = function(list, attribute) {
        if (attribute == null) {
          attribute = 'data-sort-by';
        }
        return {
          onclick: function(e) {
            var first, item;
            item = e.target.getAttribute(attribute);
            if (item) {
              first = list[0];
              list.sort(function(a, b) {
                if (isNaN(parseFloat(a[item] * 1))) {
                  if (isNaN(parseFloat(b[item] * 1))) {
                    if (a[item] > b[item]) {
                      return 1;
                    } else if (b[item] > a[item]) {
                      return -1;
                    } else {
                      return 0;
                    }
                  } else {
                    return -1;
                  }
                } else if (isNaN(parseFloat(b[item] * 1))) {
                  return 1;
                } else if (a[item] * 1 < b[item] * 1) {
                  return 1;
                } else if (b[item] * 1 < a[item] * 1) {
                  return -1;
                } else {
                  return 0;
                }
              });
              if (first === list[0]) {
                list.reverse();
              }
            }
          }
        };
      };
      return this;
    };

    Table_2.view = function(ctrl, args) {
      var begin, end, header, headers, list, offset, pageY, row, rows, x;
      rows = args.data().rows;
      headers = args.data().headers;
      pageY = state.pageY;
      begin = pageY / 60 | 0;
      end = begin + (state.pageHeight / 60 | 0 + 10);
      offset = pageY % 60;
      return m(".Occlusion", {
        style: {
          height: args.data().rows.length * 46 + "px",
          position: "relative",
          top: -offset + "px"
        }
      }, [
        m("table", {
          style: {
            top: state.pageY + "px"
          }
        }, [
          m("tr", [
            (function() {
              var k, len, ref, results1;
              ref = args.data().headers;
              results1 = [];
              for (k = 0, len = ref.length; k < len; k++) {
                header = ref[k];
                results1.push(m('th[data-sort-by=' + header + ']', ctrl.sort_table(list = args.data().rows), [header]));
              }
              return results1;
            })()
          ]), (function() {
            var k, len, ref, results1;
            ref = args.data().rows.slice(begin, +end + 1 || 9e9);
            results1 = [];
            for (x = k = 0, len = ref.length; k < len; x = ++k) {
              row = ref[x];
              results1.push(m('tr', [
                (function() {
                  var l, len1, len2, o, ref1, results2, results3;
                  if (args.checkbox) {
                    results2 = [];
                    for (l = 0, len1 = headers.length; l < len1; l++) {
                      header = headers[l];
                      if (header === "Gene") {
                        results2.push(m('td', {
                          "class": 'gene_table_item'
                        }, [
                          m('.checkbox', [
                            m('label', [
                              m('input[type=checkbox]', {
                                "class": 'checkbox gene-table-checkbox gene-search-select',
                                checked: row.selected(),
                                onclick: m.withAttr("checked", row.selected)
                              }), row[header]
                            ])
                          ])
                        ]));
                      } else {
                        results2.push(m('td', {
                          "class": 'gene_table_item'
                        }, [m('label', [row[header]])]));
                      }
                    }
                    return results2;
                  } else {
                    ref1 = args.data().headers;
                    results3 = [];
                    for (o = 0, len2 = ref1.length; o < len2; o++) {
                      header = ref1[o];
                      results3.push(m('td', [row[header]]));
                    }
                    return results3;
                  }
                })()
              ]));
            }
            return results1;
          })()
        ])
      ]);
    };

    return Table_2;

  })();

  Foo_ooo = (function(superClass) {
    extend(Foo_ooo, superClass);

    function Foo_ooo() {
      return Foo_ooo.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/foo', Foo_ooo);

    Foo_ooo.controller = function(args) {
      return this;
    };

    Foo_ooo.view = function(ctrl, args) {
      return Foo_ooo.__super__.constructor.view.call(this, m.component(Table_second));
    };

    return Foo_ooo;

  })(Page);

  Table_second = (function() {
    function Table_second() {}

    return Table_second;

  })();

  Test = (function(superClass) {
    extend(Test, superClass);

    function Test() {
      return Test.__super__.constructor.apply(this, arguments);
    }

    Routes.add("/test", Test);

    Test.controller = function(args) {
      var url;
      this.data = getEndpoint2(url = "data/meta");
      return this;
    };

    Test.view = function(ctrl) {
      return Test.__super__.constructor.view.call(this, m(".", m.component(Table_2, {
        data: ctrl.data
      })));
    };

    return Test;

  })(Page);

  dragdrop = function(element, options) {
    var activate, deactivate, update;
    activate = function(e) {
      e.preventDefault();
    };
    deactivate = function() {};
    update = function(e) {
      e.preventDefault();
      if (typeof options.onchange === 'function') {
        options.onchange((e.dataTransfer || e.target).files);
      }
    };
    options = options || {};
    element.addEventListener('dragover', activate);
    element.addEventListener('dragleave', deactivate);
    element.addEventListener('dragend', deactivate);
    element.addEventListener('drop', deactivate);
    element.addEventListener('drop', update);
    window.addEventListener('blur', deactivate);
  };

  submodule = function(module, args) {
    return module.view.bind(this, new module.controller(args));
  };

  Foo = (function() {
    function Foo() {}

    Foo.upload = function(files) {
      var formData, i, response;
      formData = new FormData;
      i = 0;
      while (i < files.length) {
        formData.append('file' + i, files[i]);
        i++;
      }
      response = m.request({
        method: 'POST',
        url: "http://" + location.hostname + ":5000/data/",
        data: formData,
        datatype: "multipart/form-data",
        serialize: function(value) {
          return value;
        },
        type: function(response) {}
      }, console.log("RESPONSE: ", JSON.stringify(response)));
      return console.log("Formdata: " + (JSON.stringify(formData.entries())));
    };

    Foo.controller = function(options) {
      options = options || {};
      this.onchange = options.onchange;
      return this;
    };

    Foo.view = function(ctrl) {
      return m('.uploader', {}, "FOO");
    };

    return Foo;

  })();

  Uploader = (function(superClass) {
    extend(Uploader, superClass);

    function Uploader() {
      return Uploader.__super__.constructor.apply(this, arguments);
    }

    Routes.add('/uploader', Uploader);

    Routes.add('/upload', Uploader);

    Uploader.controller = function() {
      this.title = function() {
        return 'Upload something';
      };
      this.foo = m.prop();
      return this;
    };

    Uploader.view = function(ctrl) {
      return m('.', m('h1', ctrl.title()), m('input[type=file]', {
        oninput: m.withAttr('value', ctrl.foo),
        value: ctrl.foo(),
        placeholder: "Username"
      }), m('input[type=button]', {
        value: "PUSH ME",
        onchange: function(files) {
          return ctrl.files = files;
        },
        onclick: function(files) {
          return alert(JSON.stringify(ctrl.foo));
        }
      }));
    };

    return Uploader;

  })(Page);


  /*
          m('.'
              m('h1', ctrl.title())
              m.component(Uploader, onchange: (files) ->
                  Uploader.upload(files)
                  #.then(->
                   *    alert 'uploaded!'
                   *    return
                  #)
                  return
              )
              m('p', 'drop a file in the box above')
          )
   */

  getEndpoint = function(url) {
    var diff, error, response;
    try {
      response = m.prop(JSON.parse(localStorage.getItem(url)));
      diff = new Date(response().date).getTime() - new Date().getTime();
      if (diff < 0) {
        throw diff;
      }
    } catch (error) {
      localStorage.setItem(url, null);
      response = m.prop(null);
    }
    if (response() === null) {
      response = m.request({
        method: "GET",
        url: "http://" + location.hostname + ":5000/" + url,
        data: {},
        datatype: 'json',
        type: function(response) {
          var binding, data, error1, item, k, l, len, len1, ref, ref1, x;
          data = {};
          data.headers = response.head.vars;
          data.rows = [];
          ref = response.results.bindings;
          for (x = k = 0, len = ref.length; k < len; x = ++k) {
            binding = ref[x];
            data.rows[x] = {};
            ref1 = response.head.vars;
            for (l = 0, len1 = ref1.length; l < len1; l++) {
              item = ref1[l];
              try {
                data.rows[x][item] = binding[item]['value'];
              } catch (error1) {
                data.rows[x][item] = '';
              }
            }
          }
          data.date = response.date;
          localStorage.setItem(url, JSON.stringify(data));
          return data;
        }
      });
    }
    return response;
  };


  /*
  This is the endpoint without front-end json manipulation.
   */

  getEndpoint2 = function(url) {
    var diff, error, response;
    try {
      response = m.prop(JSON.parse(localStorage.getItem(url)));
      diff = new Date(response().date).getTime() - new Date().getTime();
      if (diff < 0) {
        throw diff;
      }
    } catch (error) {
      localStorage.setItem(url, null);
      response = m.prop(null);
    }
    if (response() === null) {
      response = m.request({
        method: "GET",
        url: "http://" + location.hostname + ":5000/" + url,
        data: {},
        datatype: 'json',
        type: function(response) {
          localStorage.setItem(url, JSON.stringify(response));
          return response;
        }
      });
    }
    return response;
  };

  GeneSearchModel = {
    vfs: getEndpoint(url = "data/vf"),
    amrs: getEndpoint(url = "data/amr"),
    genomes: this.genomedata = getEndpoint2(url = "data/meta"),
    vfList: {},
    amrList: {},
    genomeList: [],
    selectedVF: [],
    selectedAMR: [],
    selectedGenomes: ["JHNV00000000", "ANVW00000000", "CP002729", "AMVC00000000"],
    vfresults: m.prop({}),
    amrresults: m.prop({}),
    setLists: function() {
      var gene, k, l, len, len1, ref, ref1, results1;
      this.vfList["headers"] = this.vfs().headers;
      this.amrList["headers"] = this.amrs().headers;
      ref = this.vfs().rows;
      for (k = 0, len = ref.length; k < len; k++) {
        gene = ref[k];
        gene["selected"] = m.prop(false);
        gene["visible"] = m.prop(true);
        this.vfList["rows"] = this.vfs().rows;
      }
      ref1 = this.amrs().rows;
      results1 = [];
      for (l = 0, len1 = ref1.length; l < len1; l++) {
        gene = ref1[l];
        gene["selected"] = m.prop(false);
        gene["visible"] = m.prop(true);
        results1.push(this.amrList["rows"] = this.amrs().rows);
      }
      return results1;
    },
    getSelectedVF: function() {
      var k, len, ref, results1, row;
      this.selectedVF = [];
      ref = this.vfList.rows;
      results1 = [];
      for (k = 0, len = ref.length; k < len; k++) {
        row = ref[k];
        if (row.selected()) {
          results1.push(this.selectedVF.push(row.Gene_Name));
        }
      }
      return results1;
    },
    getSelectedAMR: function() {
      var k, len, ref, results1, row;
      this.selectedAMR = [];
      ref = this.amrList.rows;
      results1 = [];
      for (k = 0, len = ref.length; k < len; k++) {
        row = ref[k];
        if (row.selected()) {
          results1.push(this.selectedAMR.push(row.Gene_Name));
        }
      }
      return results1;
    },
    reset: function() {
      var gene, k, l, len, len1, ref, ref1, results1;
      console.log("Resetting...");
      this.selectedVF = [];
      this.selectedAMR = [];
      ref = this.vfList.rows;
      for (k = 0, len = ref.length; k < len; k++) {
        gene = ref[k];
        gene["selected"] = m.prop(false);
        gene["visible"] = m.prop(true);
      }
      ref1 = this.amrList.rows;
      results1 = [];
      for (l = 0, len1 = ref1.length; l < len1; l++) {
        gene = ref1[l];
        gene["selected"] = m.prop(false);
        results1.push(gene["visible"] = m.prop(true));
      }
      return results1;
    },
    submit: function() {
      console.log("Submitting...");
      if (this.selectedVF.length === 0 && this.selectedAMR.length === 0 && this.selectedGenomes.length === 0) {
        return alert("You haven't selected any genes or genomes.");
      } else if (this.selectedGenomes.length === 0) {
        return alert("You haven't selected any genomes.");
      } else if (this.selectedVF.length === 0 && this.selectedAMR.length === 0) {
        return alert("You haven't selected any genes.");
      } else {
        this.vfresults = this.getResults(this.selectedVF);
        this.amrresults = this.getResults(this.selectedAMR);
        m.route("/results");
      }
    },
    getResults: function(geneList) {
      var response;
      response = m.prop(null);
      if (response() === null) {
        response = m.request({
          method: "POST",
          url: "http://" + location.hostname + ":5000/data/genesearchresults",
          data: {
            genome: this.selectedGenomes,
            genes: geneList
          },
          datatype: "json",
          type: function(response) {
            console.log("Results:", response);
            return response;
          }
        });
      }
      return response;
    },
    getCategories: function(type) {
      var response;
      response = m.prop(null);
      if (response() === null) {
        response = m.request({
          method: "POST",
          url: "http://" + location.hostname + ":5000/data/categories/" + type,
          data: {},
          datatype: "json",
          type: function(response) {
            console.log("Categories:", response);
            return response;
          }
        });
      }
      return response;
    }
  };

  GeneSearchModel.setLists();

  MetaData = function(args) {
    var error;
    args = args || {};
    try {
      if (this.data == null) {
        this.data = m.prop(JSON.parse(localStorage.getItem("MetaData")));
      }
    } catch (error) {
      this.data = m.prop(null);
    }
    if (this.data() === null || args.reset === true) {
      this.data = m.request({
        method: "POST",
        url: "http://" + location.hostname + ":5000/data/meta",
        data: {},
        datatype: 'json',
        type: function(response) {
          var binding, data, error1, item, k, l, len, len1, ref, ref1, x;
          data = {};
          data.headers = response.head.vars;
          data.rows = [];
          ref = response.results.bindings;
          for (x = k = 0, len = ref.length; k < len; x = ++k) {
            binding = ref[x];
            data.rows[x] = {};
            ref1 = response.head.vars;
            for (l = 0, len1 = ref1.length; l < len1; l++) {
              item = ref1[l];
              try {
                data.rows[x][item] = binding[item]['value'];
              } catch (error1) {
                data.rows[x][item] = '';
              }
            }
          }
          data.date = new Date();
          return data;
        }
      });
    }
    window.onbeforeunload = function() {
      return localStorage.setItem("MetaData", JSON.stringify(this.data()));
    };
    return this.data;
  };

  User = function(data) {
    data = data || {};
    this.username = m.prop(data.username || "");
    this.password = m.prop(data.password || "");
    this.password2 = m.prop(data.password2 || "");
    this.first_name = m.prop(data.first_name || "");
    this.last_name = m.prop(data.last_name || "");
    this.email = m.prop(data.email || "");
  };

  User.sign_up = function(user) {
    var URL, response;
    URL = 'api/users';
    return response = m.request({
      headers: get_headers(),
      method: "POST",
      url: "http://" + location.hostname + ":5000/" + URL,
      data: {
        'username': user.username(),
        'password': user.password()
      },
      datatype: 'json',
      type: function(response) {
        return response;
      }
    });
  };

  User.log_in = function(user) {
    var URL, auth, response;
    user = user || new User;
    URL = 'api/token';
    auth = "Basic " + btoa((user.username()) + ":" + (user.password()));
    return response = m.request({
      method: "GET",
      url: "http://" + location.hostname + ":5000/" + URL,
      config: function(xhr) {
        xhr.setRequestHeader('Content-Type', 'application/json');
        return xhr.setRequestHeader('Authorization', auth);
      },
      type: function(response) {
        if (typeof response === "string") {
          throw response;
        } else {
          return response;
        }
      }
    });
  };

  m.route(document.body, "/", Routes.routes);

}).call(this);
