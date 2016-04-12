class Factors extends Page
    @controller: (args) ->
        @active = m.prop("genes")
        @tabCtrl = new mc.Tabs.controller('genes')
        @vfdata = getEndpoint(url="data/vf")
        @amrdata = getEndpoint(url="data/amr")
        @genomedata = getEndpoint(url="data/meta")
        return @

    @view: (ctrl) ->
        tabOptions = {
            flavor: 'bs/nav-tabs'
            tabs: [
                { name: 'genes', label: 'Select Genes'}
                { name: 'genomes', label: 'Select Genomes'}
                { name: 'submit', label: 'Submit Query'}
            ]
        }

        renderTabContents = (ctrl) ->
            activeTab = ctrl.tabCtrl.activeTabName
            switch (activeTab())
                when 'genes' then return renderGeneForm(ctrl)
                when 'genomes' then return renderGenomeSearch(ctrl)
                when 'submit' then return renderSubmit(ctrl)
                else return m('p', "Not a valid tab name.")

        renderGeneForm = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        m.component(GeneSearchPanel, {
                            title: "Virulence Factor"
                            type: "vf"
                            data: ctrl.genedata
                        })
                        m.component(GeneSearchPanel, {
                            title: "Antimicrobial Resistance"
                            type: "amr"
                            data: ctrl.genedata
                        })
                    ])
                ])
            ])

        renderGenomeSearch = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        "Testing Genome Search"
                    ])
                ])
            ])

        renderSubmit = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        "Testing Submit"
                    ])
                ])
            ])


        super(
            m('#wrapper', [
                m.component(Sidebar)
                m('.', {id: 'page-content-wrapper'}, [
                    m('.', {id: 'page-content -inset'}, [
                        m('.container-fluid', [
                            m.component(FactorsIntro)
                            m('.container', [
                                mc.Tabs.view(ctrl.tabCtrl, tabOptions)
                                renderTabContents(ctrl)
                            ])
                        ])
                    ])
                ])
            ])
        )

FactorsIntro =
    view: (ctrl, args) ->
        m(".intro", [
            m('.', {class: 'row'}, [
                m('.col-xs-8', [
                    m('.content-header', [
                        m('h1', [
                            m('span', {class: "title_part1"}, 'VIRULENCE & AMR '),
                            m('span', {class: "title_part2"}, 'GENES')
                        ])
                    ])
                ])
                m('.', {class: 'col-xs-4'}, [
                    m('button', {id: "intro-button", \
                                 class: "btn btn-danger btn-lg", \
                                 type:"button"}, [
                        "INTRODUCTION"
                    ])
                ])
            ])
            m("p", "Search for the presence or absence of
                    virulence factor genes or antimicrobial
                    resistance genes in genomes of interest.
                    Detailed information on individual virulence
                    factor or antimicrobial resistance genes
                    can be retrieved by clicking on the
                    individual genes.")
        ])

## Model that holds the gene and genome data
##  and the selected factors/genomes.
DataModel =
    vfList: []
    amrList: []
    genomeslist: []
    selectedVF: []
    selectedAMR: []
    selectedGenomes: []


