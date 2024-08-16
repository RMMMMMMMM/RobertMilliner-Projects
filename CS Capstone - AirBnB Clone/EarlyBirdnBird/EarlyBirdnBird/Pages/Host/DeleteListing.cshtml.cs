using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;

namespace EarlyBirdnBird.Pages.Host
{
    public class DeleteListingModel : PageModel
    {


        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;


        //DB injection
        private readonly UnitOfWork _uow;
        //second injection to make the file path possible for images know where www.root is located
        private readonly IWebHostEnvironment _webHostEnvironment;



        public int rentalID { get; set; }
        public string rentalNAME { get; set; }


        public DeleteListingModel(ApplicationDbContext db, UserManager<AppUser> um, UnitOfWork uow, IWebHostEnvironment webHostEnvironment)
        {
            _dbContext = db;
            _userManager = um;
            _uow = uow;
            _webHostEnvironment = webHostEnvironment;
        }





        public IActionResult OnGetAsync()
        {

            var temp = Request.Query["rentalID"];

            rentalID = Convert.ToInt32(temp);


            var tempRemtal = (from r in _dbContext.Rentals
                              where r.Id == rentalID
                              select r).FirstOrDefault();

            rentalNAME = tempRemtal.Name;


            return Page();
        }


        public async Task<IActionResult> OnPostAsync()
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }

            var temp = Request.Query["rentalID"];

            rentalID = Convert.ToInt32(temp);


            var tempRemtal = (from r in _dbContext.Rentals
                              where r.Id == rentalID
                              select r).FirstOrDefault();

            //if(tempRemtal.Images != null && tempRemtal.Images.Any())
            //{
            //    var image = await _dbContext.Images.FirstOrDefaultAsync(i => i.RentalId == tempRemtal.Id);
            //    _dbContext.Images.Remove(image);
            //    await _dbContext.SaveChangesAsync();}



            

                // Find all associated image records
                var images = await _dbContext.Images.Where(i => i.RentalId == rentalID).ToListAsync();

                // Remove all image records from the context
                _dbContext.Images.RemoveRange(images);

                // Remove the rental record from the context
                _dbContext.Rentals.Remove(tempRemtal);
                   
            await _dbContext.SaveChangesAsync();


            return RedirectToPage("./Listings");


        }


    }
}
