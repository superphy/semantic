###
Mock_component.coffee
    
    A trial table to test having a static table view and multiple data views in one page
    URL: localhost/superphy/?/?????

Args passed in:
    debug_data: Data to populate the table
###
# coffeelint: disable=max_line_length

###
    MOCK_COMPONENT:
###
class Mockcomponent
    
    TableView = (TableData) ->
        this.row = TableData.row or m.prop([]) # The m.prop() here is an error handler, if 
        this.headers = TableData.headers

        m("table", {"border": "2px solid black", "table-layout": "fixed", "overflow": "auto"}
            m("tr", [
                for header in headers
                    m('th', {width: "100px"}, header)
            ])
            
            for row, x in this.row #when JSON.stringify(row).search(/Unknown/i) > -1
                m('tr' 
                    for column in this.headers
                        m('td', {width: "100px"}, [row[column]])
                )
        )
 

    ###
        DATAMOCK:
        nested class for getting static mockdata
    ###
    class DataMock
        constructor: ->
            this.row = [{"Accession":"AAJX00000000","Bioproject_Id":"15630","Biosample_Id":"2435887","Common_Name":"","Genome_Uri":"https://github.com/superphy#AAJX00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"B171","Syndromes":""},{"Accession":"ADTX00000000","Bioproject_Id":"47257","Biosample_Id":"189171","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ADTX00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"BEI RESOURES STRAIN HM-351,\nMS 153-1","Syndromes":""},{"Accession":"ADUO00000000","Bioproject_Id":"40269","Biosample_Id":"792452","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ADUO00000000","Geographic_Location":"Bangladesh","Isolation_Date":"1983-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"2","Serotype_O":"114","Strain":"AEPEC,\nE128010","Syndromes":""},{"Accession":"AERP00000000","Bioproject_Id":"61463","Biosample_Id":"2472077","Common_Name":"","Genome_Uri":"https://github.com/superphy#AERP00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"7","Serotype_O":"157","Strain":"A,\n1044","Syndromes":""},{"Accession":"AEZK00000000","Bioproject_Id":"51089","Biosample_Id":"2436554","Common_Name":"cow","Genome_Uri":"https://github.com/superphy#AEZK00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"Bos taurus","Serotype_H":"","Serotype_O":"","Strain":"5.0588","Syndromes":""},{"Accession":"AIBQ00000000","Bioproject_Id":"79303","Biosample_Id":"2435993","Common_Name":"human","Genome_Uri":"https://github.com/superphy#AIBQ00000000","Geographic_Location":"Slovenia","Isolation_Date":"","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"C79_08","Syndromes":""},{"Accession":"AIGR00000000","Bioproject_Id":"51003","Biosample_Id":"190940","Common_Name":"","Genome_Uri":"https://github.com/superphy#AIGR00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"DEC10C","Syndromes":""},{"Accession":"AIHG00000000","Bioproject_Id":"51033","Biosample_Id":"199318","Common_Name":"","Genome_Uri":"https://github.com/superphy#AIHG00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"DEC13B","Syndromes":""},{"Accession":"AIHM00000000","Bioproject_Id":"51047","Biosample_Id":"210648","Common_Name":"","Genome_Uri":"https://github.com/superphy#AIHM00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"DEC14C","Syndromes":""},{"Accession":"AJVR00000000","Bioproject_Id":"129421","Biosample_Id":"991747","Common_Name":"","Genome_Uri":"https://github.com/superphy#AJVR00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"2","Serotype_O":"103","Strain":"CVM9450","Syndromes":""},{"Accession":"AKMJ00000000","Bioproject_Id":"65931","Biosample_Id":"2435930","Common_Name":"human","Genome_Uri":"https://github.com/superphy#AKMJ00000000","Geographic_Location":"","Isolation_Date":"2006-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"7","Serotype_O":"157","Strain":"EC4439","Syndromes":""},{"Accession":"AMTO00000000","Bioproject_Id":"65751","Biosample_Id":"1036883","Common_Name":"","Genome_Uri":"https://github.com/superphy#AMTO00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"10.0821","Syndromes":""},{"Accession":"AMVC00000000","Bioproject_Id":"65969","Biosample_Id":"808713","Common_Name":"","Genome_Uri":"https://github.com/superphy#AMVC00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"EC1865","Syndromes":""},{"Accession":"ANTT00000000","Bioproject_Id":"163461","Biosample_Id":"974079","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ANTT00000000","Geographic_Location":"","Isolation_Date":"2010-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"KTE224","Syndromes":""},{"Accession":"ANVW00000000","Bioproject_Id":"164875","Biosample_Id":"854618","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ANVW00000000","Geographic_Location":"Denmark","Isolation_Date":"2010-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"KTE59","Syndromes":""},{"Accession":"ANWG00000000","Bioproject_Id":"164967","Biosample_Id":"854664","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ANWG00000000","Geographic_Location":"Denmark","Isolation_Date":"2010-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"KTE123","Syndromes":""},{"Accession":"AOEM00000000","Bioproject_Id":"65893","Biosample_Id":"1112702","Common_Name":"","Genome_Uri":"https://github.com/superphy#AOEM00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"PA47","Syndromes":""},{"Accession":"","Bioproject_Id":"","Biosample_Id":"","Common_Name":"","Genome_Uri":"https://github.com/superphy#AP009048","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"","Syndromes":""},{"Accession":"APWX00000000","Bioproject_Id":"183614","Biosample_Id":"2732275","Common_Name":"human","Genome_Uri":"https://github.com/superphy#APWX00000000","Geographic_Location":"Germany","Isolation_Date":"2000-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"21","Serotype_O":"174","Strain":"220/00","Syndromes":""},{"Accession":"APWZ00000000","Bioproject_Id":"77447","Biosample_Id":"829231","Common_Name":"","Genome_Uri":"https://github.com/superphy#APWZ00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"178850","Syndromes":""},{"Accession":"ATOC00000000","Bioproject_Id":"203078","Biosample_Id":"2436903","Common_Name":"human","Genome_Uri":"https://github.com/superphy#ATOC00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"AB42554418-ISOLATE1","Syndromes":"Bacteremia"},{"Accession":"AVRM00000000","Bioproject_Id":"72911","Biosample_Id":"1919630","Common_Name":"","Genome_Uri":"https://github.com/superphy#AVRM00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"T1282_01","Syndromes":""},{"Accession":"AVWR00000000","Bioproject_Id":"186153","Biosample_Id":"1885777","Common_Name":"human","Genome_Uri":"https://github.com/superphy#AVWR00000000","Geographic_Location":"Denmark","Isolation_Date":"2004-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"HVH 132 4-6876862","Syndromes":"Bacteremia,\nUrinary tract infection (cystitis)"},{"Accession":"AWAB00000000","Bioproject_Id":"186257","Biosample_Id":"1885867","Common_Name":"human","Genome_Uri":"https://github.com/superphy#AWAB00000000","Geographic_Location":"Denmark","Isolation_Date":"2005-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"KOEGE 10 25A","Syndromes":"Bacteremia,\nUrinary tract infection (cystitis)"},{"Accession":"AWQN00000000","Bioproject_Id":"215830","Biosample_Id":"2335231","Common_Name":"","Genome_Uri":"https://github.com/superphy#AWQN00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"WA1","Syndromes":""},{"Accession":"AWQO00000000","Bioproject_Id":"215830","Biosample_Id":"2335232","Common_Name":"","Genome_Uri":"https://github.com/superphy#AWQO00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"WA2","Syndromes":""},{"Accession":"AYHV00000000","Bioproject_Id":"202049","Biosample_Id":"2138668","Common_Name":"","Genome_Uri":"https://github.com/superphy#AYHV00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"BIDMC 37","Syndromes":""},{"Accession":"","Bioproject_Id":"","Biosample_Id":"","Common_Name":"","Genome_Uri":"https://github.com/superphy#BA000007","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"","Syndromes":""},{"Accession":"CCQB00000000","Bioproject_Id":"260147","Biosample_Id":"3090438","Common_Name":"","Genome_Uri":"https://github.com/superphy#CCQB00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"","Serotype_O":"","Strain":"PUERTO RICAN","Syndromes":""},{"Accession":"CP000247","Bioproject_Id":"16235","Biosample_Id":"2604181","Common_Name":"","Genome_Uri":"https://github.com/superphy#CP000247","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"31","Serotype_O":"6","Strain":"536","Syndromes":""},{"Accession":"CP002729","Bioproject_Id":"42137","Biosample_Id":"2604221","Common_Name":"pig","Genome_Uri":"https://github.com/superphy#CP002729","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"Sus scrofa","Serotype_H":"","Serotype_O":"149","Strain":"UMNK88","Syndromes":"Diarrhea"},{"Accession":"CP003109","Bioproject_Id":"68245","Biosample_Id":"2604278","Common_Name":"human","Genome_Uri":"https://github.com/superphy#CP003109","Geographic_Location":"United States, California","Isolation_Date":"1974-12-01","Scientific_Name":"Homo sapiens","Serotype_H":"7","Serotype_O":"55","Strain":"RM12579","Syndromes":""},{"Accession":"JHFM00000000","Bioproject_Id":"218110","Biosample_Id":"2352975","Common_Name":"","Genome_Uri":"https://github.com/superphy#JHFM00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"19","Serotype_O":"121","Strain":"2010C-3609","Syndromes":""},{"Accession":"JHJJ00000000","Bioproject_Id":"218110","Biosample_Id":"2353082","Common_Name":"","Genome_Uri":"https://github.com/superphy#JHJJ00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"41","Serotype_O":"169","Strain":"F9792","Syndromes":""},{"Accession":"JHNI00000000","Bioproject_Id":"218110","Biosample_Id":"2352910","Common_Name":"","Genome_Uri":"https://github.com/superphy#JHNI00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"7","Serotype_O":"157","Strain":"06-3745","Syndromes":""},{"Accession":"JHNV00000000","Bioproject_Id":"218110","Biosample_Id":"2352897","Common_Name":"","Genome_Uri":"https://github.com/superphy#JHNV00000000","Geographic_Location":"","Isolation_Date":"","Scientific_Name":"","Serotype_H":"4","Serotype_O":"119","Strain":"03-3458","Syndromes":""},{"Accession":"JMIZ00000000","Bioproject_Id":"245436","Biosample_Id":"2739483","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JMIZ00000000","Geographic_Location":"United States, California, Davis","Isolation_Date":"2013-01-01","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"UCD_JA03_PB","Syndromes":""},{"Accession":"JNMX00000000","Bioproject_Id":"233910","Biosample_Id":"2680260","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JNMX00000000","Geographic_Location":"Tanzania","Isolation_Date":"2009-07-08","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"3-475-03_S4_C2","Syndromes":""},{"Accession":"JNNK00000000","Bioproject_Id":"246239","Biosample_Id":"2744160","Common_Name":"","Genome_Uri":"https://github.com/superphy#JNNK00000000","Geographic_Location":"Germany","Isolation_Date":"1992-01-01","Scientific_Name":"","Serotype_H":"49","Serotype_O":"181","Strain":"S.93","Syndromes":""},{"Accession":"JNPI00000000","Bioproject_Id":"233778","Biosample_Id":"2687494","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JNPI00000000","Geographic_Location":"Tanzania","Isolation_Date":"2009-04-24","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"2-177-06_S3_C1","Syndromes":""},{"Accession":"JNQR00000000","Bioproject_Id":"233839","Biosample_Id":"2680240","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JNQR00000000","Geographic_Location":"Tanzania","Isolation_Date":"2009-04-24","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"2-316-03_S3_C3","Syndromes":""},{"Accession":"JNRT00000000","Bioproject_Id":"233933","Biosample_Id":"2680205","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JNRT00000000","Geographic_Location":"Tanzania","Isolation_Date":"2009-07-08","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"3-020-07_S4_C3","Syndromes":""},{"Accession":"JOSM00000000","Bioproject_Id":"233780","Biosample_Id":"2689034","Common_Name":"human","Genome_Uri":"https://github.com/superphy#JOSM00000000","Geographic_Location":"Tanzania","Isolation_Date":"2009-04-23","Scientific_Name":"Homo sapiens","Serotype_H":"","Serotype_O":"","Strain":"2-222-05_S3_C1","Syndromes":""}]
            this.headers = ["Genome_Uri","Syndromes","Accession","Biosample_Id","Bioproject_Id","Strain","Serotype_O","Serotype_H","Scientific_Name","Common_Name","Isolation_Date","Geographic_Location"]
            this.len = this.row.length

    ###
        PageFunction: 
        Pagination class for modelling the change of pages
    ###
    class PageFunction
        ###
        constructor: (offs, len)->
            this.offs = m.prop(offs) # This is the number of row displayed at a time
            this.currentp = m.prop(1) # To display next to next and previous buttons
            this.numRows = m.prop(this.offs()-1) # A Class because ofm.prop. Sets the upper limit for number of rows to display at the time
            this.start = m.prop(0) # Row to begin display, adjusted in next()
            this.pageNumber = Math.floor(len/this.offs())+1 # This divides the number of records by the number of rows dispayed, and gets the number of pages
            console.log("This is the number of pages #{this.pageNumber}")
        nextpage: (max)->
            if (max == false)
                return # The last row has been printed, do nothing
            else
                this.start(this.start() + this.offs()) # Started at 0, increments by 50
                this.currentp(this.currentp() + 1)
                this.numRows(this.numRows() + this.offs())
        prevpage: ->
            if (this.currentp() == 1)
               return
            else
                this.numRows(this.numRows() - this.offs()) # Started at 50, decrements by 50 to change the offs in data stream
                this.start(this.start() - this.offs()) # Started at 0, decrements by 50 
                this.currentp(this.currentp() - 1)
        setpage: (num) ->
            this.currentp = num
        createButtons: (mockdata) ->
            return
        ###


    # Mock_component controller
    @controller: () ->
        this.mockdata = new DataMock
        this.pagebuttons = new PageFunction(25, this.mockdata.len)
        console.log("This is the number of rows #{this.mockdata.len}")
        return @

    @view: (ctrl) -> 
        m('.'
            #m('.mapdiv', {config: loadmap})
            #m('map', {config: loadmap})
            TableView(ctrl.mockdata),
            #ctrl.pagebuttons.createButtons(ctrl.mockdata)
        )








    ### Tested out the map view, not working yet
    loadmap = (map, init) ->
        if init
            return
        else
            map = new google.maps.Map(document.getElementById('map'), {
                center: {lat: 49.69349, lng: -112.84183999999999},
                zoom: 4
            });
        
        #map = new (google.maps.Map)(element, opts)
        #opts = {
        #    center: new (google.maps.LatLng)(49.69349,-112.84183999999999)
        #    zoom: 4
        #    mapTypeId: google.maps.MapTypeId.ROADMAP
        #}
    ###

    ###
    m('.'
        TableView(ctrl.mockdata),
        #pagebuttons.view,
        m("button.btn btn-primary[type=button]" 
            #onclick: -> ctrl.pagebuttons.setpage(m.withattr(value)) }
            #value: ctrl.pagebuttons.setpage()
            "Next Page"
        )
    )
    #return m.component(Table_2, ctrl.mockdata)
    ###


