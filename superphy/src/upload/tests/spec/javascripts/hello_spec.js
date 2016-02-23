/**
 * Created by Stephen Kan on 20/10/15.
 Example for making javascript tests
 */

describe("Test hello.js", function() {
    it("return 'hello'", function () {
        var result = "Hello"//Hello.world();
        expect(result).toBe("Hello");
    });
});