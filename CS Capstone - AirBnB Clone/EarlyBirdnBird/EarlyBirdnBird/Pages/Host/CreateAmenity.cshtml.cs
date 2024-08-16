using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using System.ComponentModel.DataAnnotations;


//This page is not currently used
//keeping because we could use it 
//maybe as an Admin feature


namespace EarlyBirdnBird.Pages.Host
{
    public class CreateAmenityModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;
        


        public CreateAmenityModel(ApplicationDbContext db, UserManager<AppUser> um)
        {
            _dbContext = db;
            _userManager = um;
        }


        [BindProperty]
        public InputModel Input { get; set; }

        public class InputModel
        {
            [Required]
            [DataType(DataType.Text)]
            [Display(Name = "Name")]
            public string Name { get; set; }

        }



        public IActionResult OnGet()
        {
            return Page();
        }


        public async Task<IActionResult> OnPostAsync()
        {

            Amenity newAmenity = new Amenity()
            {
                Name = Input.Name
            };
            _dbContext.Amenities.Add(newAmenity);
            await _dbContext.SaveChangesAsync();

            //placeholder return so program builds/runs
            return RedirectToPage("./Listings");
        }


    }
}
