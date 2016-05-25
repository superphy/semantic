###
Table_Second.coffee
    A trial table to test a different pagination style.

Args passed in:
    data: Data to populate the table
###

# coffeelint: disable=max_line_length

class Table_second
# Class to define a table with standard pagination (next and previous arrows), and optionally dropdown to select the number of rows to be displayed
    state = {pageY: 0, pageHeight: window.innerHeight}
    window.addEventListener("scroll", (e) ->
        state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
        state.pageHeight = window.innerHeight
        m.redraw()
    ) # This is declaring the window height and width

    Routes.add('/tablesecond', @) # The URL of the page localhost/superphy/?/tablesecond
    
    @controller: (args) -> 
        ###
        vm.add = function() {
            if (vm.description()) {
                vm.list.push(new todo.Todo({description: vm.description()}));
                vm.description("");
            }
        };
        ###
    
    
        this.userChoice = prompt('Do you choose rock, paper, or scissors?')
        this.computerChoice = Math.random()
        this.result = 'result'
        
        # To check for a tie
        if (this.computerChoice < 0.34)
            this.computerChoice = 'rock'
        else if (this.computerChoice <= 0.67)
            this.computerChoice = 'paper'
        else
            this.computerChoice = 'scissors'
            console.log('You chose: ' + this.userChoice)
            console.log('Computer chose: ' + this.computerChoice + '\n')
        # IF IT A TIE
        
        @compare = (choice1, choice2) ->

            console.log("Choice1 = " + choice1)
            console.log("Choice2  = " + choice2)

            if(choice1 == choice2)
                this.userChoice = prompt('That was a tie! \n\nPlease choose again :)')
                if this.computerChoice < 0.34
                    this.computerChoice = 'rock'
                else if this.computerChoice <= 0.67
                    this.computerChoice = 'paper'
                else
                    this.computerChoice = 'scissors'
    
                console.log('_________________________________')
                console.log("Your NEW choice: #{this.userChoice}")
                console.log('Computer NEW choice: ' + this.computerChoice + ' \n')
                this.compare(this.userChoice, this.computerChoice)

            # IF CHOICE IS ROCK
            else if(choice1 == 'rock') 
                if choice2 == 'scissors'
                    console.log 'rock wins'
                    this.result = 'rock wins'
                else
                    console.log 'paper wins'
            # IF CHOICE IS PAPER
            else if(choice1 == 'paper')
                if choice2 == 'rock'
                    console.log 'paper wins'
                    this.result = 'paper wins'
                else
                    console.log 'scissors wins'
                    this.result = 'scissors wins'
            # IF CHOICE IS SCISSORS
            else if(choice1 == 'scissors')
                if choice2 == 'paper'
                    console.log 'scissors wins'
                    this.result = 'scissors wins'
            else
                  console.log 'rock wins'
                  this.result = 'rock wins'
            # end else-if scissors
            # return
            # End compare function 

        return @

    @view: (ctrl, args) ->
        
        res = ctrl.result
        ctrl.compare(ctrl.userChoice, ctrl.computerChoice) 
        
        #m('h1', ['The Winner is: ', ctrl.result, '!'])
        #m('h1', {style: 'color:blue'}, "TABLE SECOND!!! PLEASE PRINT!! " + ctrl.res)
        #rows = args.data().rows
        #headers = args.data().headers

        pageY = state.pageY
        begin = pageY / 60 | 0
        
        return m("html", [ #first bracket
            m("body", [
                m('h1', ['The winner? ', ctrl.result, '!']),
                m("input"),
                m("button", {onclick: -> alert('FOO')}, "Create"),
                m("table", {style: {top: state.pageY + "px"}}, [
                    m("tr", [
                        m("td", [
                            m("input[type=checkbox]")
                        ]),
                        m("td", "task description"),
                    ])
                ])
            ])
        ]); #Closing 

        # TUESDAY: Investigate what the "html" tag means
        #   - TRY this in the Table view code
        #   - DATA flow, like the /test page